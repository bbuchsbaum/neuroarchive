#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

const int64_t Q_MIN = -2147483648LL;
const int64_t Q_MAX = 2147483647LL;
const uint8_t WIDTH_MAX = 34;

uint64_t checked_limit(double value, const char* name) {
  if (!R_finite(value) || value < 1.0 || value != std::floor(value) ||
      value > 9007199254740991.0) {
    stop(std::string(name) + " must be a positive exact integer");
  }
  return static_cast<uint64_t>(value);
}

uint64_t checked_product(uint64_t a, uint64_t b, const char* what) {
  if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a) {
    stop(std::string(what) + " overflows uint64");
  }
  return a * b;
}

int64_t checked_q(double value) {
  if (!R_finite(value) || value != std::floor(value) ||
      value < static_cast<double>(Q_MIN) || value > static_cast<double>(Q_MAX)) {
    stop("q must contain finite whole numbers in the signed-int32 domain");
  }
  return static_cast<int64_t>(value);
}

uint64_t zigzag(int64_t value) {
  if (value >= 0) return static_cast<uint64_t>(value) << 1;
  return (static_cast<uint64_t>(-(value + 1)) << 1) | 1ULL;
}

int64_t unzigzag(uint64_t value) {
  if ((value & 1ULL) == 0) return static_cast<int64_t>(value >> 1);
  return -static_cast<int64_t>(value >> 1) - 1;
}

uint8_t bit_width(const std::vector<uint64_t>& values) {
  uint64_t maximum = 0;
  for (size_t i = 0; i < values.size(); ++i) maximum = std::max(maximum, values[i]);
  uint8_t width = 0;
  while (maximum != 0) {
    ++width;
    maximum >>= 1;
  }
  return width;
}

struct PredictorPlan {
  uint8_t id;
  std::vector<int64_t> anchors;
  std::vector<uint64_t> unsigned_residuals;
  uint8_t width;
  uint64_t cost_bits;
};

PredictorPlan make_plan(const std::vector<int64_t>& q, uint8_t id) {
  PredictorPlan plan;
  plan.id = id;
  if (id == 0) {
    for (size_t t = 0; t < q.size(); ++t) plan.unsigned_residuals.push_back(zigzag(q[t]));
  } else if (id == 1) {
    plan.anchors.push_back(q[0]);
    for (size_t t = 1; t < q.size(); ++t) {
      plan.unsigned_residuals.push_back(zigzag(q[t] - q[t - 1]));
    }
  } else {
    plan.anchors.push_back(q[0]);
    plan.anchors.push_back(q[1]);
    for (size_t t = 2; t < q.size(); ++t) {
      const int64_t residual = q[t] - 2LL * q[t - 1] + q[t - 2];
      plan.unsigned_residuals.push_back(zigzag(residual));
    }
  }
  plan.width = bit_width(plan.unsigned_residuals);
  if (plan.width > WIDTH_MAX) stop("internal QPP residual width exceeds 34 bits");
  plan.cost_bits = static_cast<uint64_t>(plan.anchors.size()) * 32ULL +
    static_cast<uint64_t>(plan.unsigned_residuals.size()) * plan.width;
  return plan;
}

PredictorPlan choose_plan(const std::vector<int64_t>& q) {
  PredictorPlan best = make_plan(q, 0);
  PredictorPlan p1 = make_plan(q, 1);
  if (p1.cost_bits < best.cost_bits) best = p1;
  if (q.size() >= 2) {
    PredictorPlan p2 = make_plan(q, 2);
    if (p2.cost_bits < best.cost_bits) best = p2;
  }
  return best;
}

class BoundedBytes {
 public:
  explicit BoundedBytes(uint64_t limit) : limit_(limit) {}

  void byte(uint8_t value) {
    ensure(1);
    data_.push_back(value);
  }

  void u16(uint16_t value) {
    ensure(2);
    data_.push_back(static_cast<uint8_t>(value & 0xff));
    data_.push_back(static_cast<uint8_t>((value >> 8) & 0xff));
  }

  void u32(uint32_t value) {
    ensure(4);
    for (int shift = 0; shift < 32; shift += 8) {
      data_.push_back(static_cast<uint8_t>((value >> shift) & 0xff));
    }
  }

  void i32(int64_t value) {
    if (value < Q_MIN || value > Q_MAX) stop("anchor is outside signed-int32 domain");
    const uint32_t raw = value < 0
      ? static_cast<uint32_t>(static_cast<uint64_t>(value + 4294967296LL))
      : static_cast<uint32_t>(value);
    u32(raw);
  }

  void append(const std::vector<uint8_t>& values) {
    ensure(values.size());
    data_.insert(data_.end(), values.begin(), values.end());
  }

  const std::vector<uint8_t>& data() const { return data_; }

 private:
  void ensure(uint64_t extra) {
    if (extra > limit_ || static_cast<uint64_t>(data_.size()) > limit_ - extra) {
      stop("encoded QPP tile exceeds max_bytes");
    }
  }
  uint64_t limit_;
  std::vector<uint8_t> data_;
};

std::vector<uint8_t> pack_lsb(const std::vector<uint64_t>& values, uint8_t width) {
  if (width == 0) return std::vector<uint8_t>();
  const uint64_t bits = checked_product(values.size(), width, "packed bit count");
  const uint64_t nbyte = (bits + 7ULL) / 8ULL;
  if (nbyte > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
    stop("packed payload is too large for this platform");
  }
  std::vector<uint8_t> out(static_cast<size_t>(nbyte), 0);
  uint64_t cursor = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    const uint64_t value = values[i];
    if (width < 64 && value >= (1ULL << width)) stop("residual does not fit width");
    for (uint8_t bit = 0; bit < width; ++bit, ++cursor) {
      if ((value >> bit) & 1ULL) {
        out[static_cast<size_t>(cursor >> 3)] |= static_cast<uint8_t>(1U << (cursor & 7ULL));
      }
    }
  }
  return out;
}

class ByteReader {
 public:
  explicit ByteReader(const RawVector& bytes) : bytes_(bytes), pos_(0) {}

  uint8_t byte() {
    require(1);
    return bytes_[pos_++];
  }

  uint16_t u16() {
    require(2);
    uint16_t value = static_cast<uint16_t>(bytes_[pos_]) |
      static_cast<uint16_t>(bytes_[pos_ + 1]) << 8;
    pos_ += 2;
    return value;
  }

  uint32_t u32() {
    require(4);
    uint32_t value = 0;
    for (int shift = 0; shift < 32; shift += 8) {
      value |= static_cast<uint32_t>(bytes_[pos_++]) << shift;
    }
    return value;
  }

  int64_t i32() {
    const uint32_t value = u32();
    if (value <= 2147483647U) return static_cast<int64_t>(value);
    return static_cast<int64_t>(value) - 4294967296LL;
  }

  size_t position() const { return pos_; }
  size_t size() const { return static_cast<size_t>(bytes_.size()); }

  const RawVector& bytes() const { return bytes_; }

  void advance(size_t n) {
    require(n);
    pos_ += n;
  }

 private:
  void require(size_t n) const {
    if (n > size() || pos_ > size() - n) stop("truncated QPP tile");
  }
  const RawVector& bytes_;
  size_t pos_;
};

std::vector<uint64_t> unpack_lsb(ByteReader& reader, uint32_t count,
                                 uint8_t width, uint32_t payload_nbyte) {
  const uint64_t used_bits = checked_product(count, width, "residual bit count");
  const uint64_t expected_nbyte = (used_bits + 7ULL) / 8ULL;
  if (expected_nbyte != payload_nbyte) stop("inconsistent QPP packed payload length");
  const size_t start = reader.position();
  if (payload_nbyte > reader.size() || start > reader.size() - payload_nbyte) {
    stop("truncated QPP tile");
  }
  const RawVector& bytes = reader.bytes();
  std::vector<uint64_t> values(count, 0);
  uint64_t cursor = 0;
  for (uint32_t i = 0; i < count; ++i) {
    uint64_t value = 0;
    for (uint8_t bit = 0; bit < width; ++bit, ++cursor) {
      const uint8_t packed = bytes[start + static_cast<size_t>(cursor >> 3)];
      value |= static_cast<uint64_t>((packed >> (cursor & 7ULL)) & 1U) << bit;
    }
    values[i] = value;
  }
  for (uint64_t bit = used_bits; bit < static_cast<uint64_t>(payload_nbyte) * 8ULL; ++bit) {
    const uint8_t packed = bytes[start + static_cast<size_t>(bit >> 3)];
    if (((packed >> (bit & 7ULL)) & 1U) != 0) stop("non-zero QPP padding bits");
  }
  reader.advance(payload_nbyte);
  return values;
}

void assign_checked(NumericMatrix& out, int row, int col, int64_t value) {
  if (value < Q_MIN || value > Q_MAX) stop("decoded q overflows signed-int32 domain");
  out(row, col) = static_cast<double>(value);
}

} // namespace

// [[Rcpp::export]]
RawVector qpp_encode_tile_rcpp(NumericMatrix q, double max_values,
                               double max_bytes) {
  const uint64_t values_limit = checked_limit(max_values, "max_values");
  const uint64_t bytes_limit = checked_limit(max_bytes, "max_bytes");
  const int n_time = q.nrow();
  const int n_parcel = q.ncol();
  if (n_time < 1 || n_time > 65535 || n_parcel < 1 || n_parcel > 65535) {
    stop("QPP tile dimensions must be in 1..65535");
  }
  const uint64_t value_count = checked_product(n_time, n_parcel, "QPP value count");
  if (value_count > values_limit) stop("QPP tile exceeds max_values");

  BoundedBytes out(bytes_limit);
  out.byte('Q'); out.byte('P'); out.byte('P'); out.byte('1');
  out.u16(static_cast<uint16_t>(n_time));
  out.u16(static_cast<uint16_t>(n_parcel));

  std::vector<int64_t> column(static_cast<size_t>(n_time));
  for (int k = 0; k < n_parcel; ++k) {
    for (int t = 0; t < n_time; ++t) column[t] = checked_q(q(t, k));
    const PredictorPlan plan = choose_plan(column);
    const std::vector<uint8_t> payload = pack_lsb(plan.unsigned_residuals, plan.width);
    if (payload.size() > std::numeric_limits<uint32_t>::max()) {
      stop("QPP parcel payload exceeds uint32 length");
    }
    out.byte(plan.id);
    out.byte(plan.width);
    out.byte(static_cast<uint8_t>(plan.anchors.size()));
    out.byte(0);
    out.u16(static_cast<uint16_t>(plan.unsigned_residuals.size()));
    out.u32(static_cast<uint32_t>(payload.size()));
    for (size_t i = 0; i < plan.anchors.size(); ++i) out.i32(plan.anchors[i]);
    out.append(payload);
  }
  RawVector result(out.data().size());
  std::copy(out.data().begin(), out.data().end(), result.begin());
  return result;
}

// [[Rcpp::export]]
NumericMatrix qpp_decode_tile_rcpp(RawVector bytes, int max_time,
                                   int max_parcels, double max_values,
                                   double max_bytes) {
  if (max_time < 1 || max_parcels < 1) stop("decode dimension limits must be positive");
  const uint64_t values_limit = checked_limit(max_values, "max_values");
  const uint64_t bytes_limit = checked_limit(max_bytes, "max_bytes");
  if (static_cast<uint64_t>(bytes.size()) > bytes_limit) stop("QPP tile exceeds max_bytes");
  ByteReader reader(bytes);
  if (reader.byte() != 'Q' || reader.byte() != 'P' ||
      reader.byte() != 'P' || reader.byte() != '1') stop("invalid QPP magic");
  const uint16_t n_time = reader.u16();
  const uint16_t n_parcel = reader.u16();
  if (n_time < 1 || n_parcel < 1) stop("impossible QPP dimensions");
  if (n_time > max_time || n_parcel > max_parcels) stop("QPP dimensions exceed decode limits");
  const uint64_t value_count = checked_product(n_time, n_parcel, "QPP value count");
  if (value_count > values_limit) stop("QPP tile exceeds max_values");

  NumericMatrix q(n_time, n_parcel);
  for (int k = 0; k < n_parcel; ++k) {
    const uint8_t predictor = reader.byte();
    const uint8_t width = reader.byte();
    const uint8_t anchor_count = reader.byte();
    const uint8_t reserved = reader.byte();
    const uint16_t residual_count = reader.u16();
    const uint32_t payload_nbyte = reader.u32();
    if (predictor > 2 || width > WIDTH_MAX || reserved != 0) stop("invalid QPP header");
    const uint8_t expected_anchor = predictor;
    if (anchor_count != expected_anchor || anchor_count > n_time ||
        residual_count != static_cast<uint16_t>(n_time - anchor_count)) {
      stop("inconsistent QPP predictor lengths");
    }
    std::vector<int64_t> anchors(anchor_count);
    for (uint8_t i = 0; i < anchor_count; ++i) anchors[i] = reader.i32();
    const std::vector<uint64_t> packed = unpack_lsb(
      reader, residual_count, width, payload_nbyte
    );
    std::vector<int64_t> residuals(packed.size());
    for (size_t i = 0; i < packed.size(); ++i) residuals[i] = unzigzag(packed[i]);

    if (predictor == 0) {
      for (int t = 0; t < n_time; ++t) assign_checked(q, t, k, residuals[t]);
    } else if (predictor == 1) {
      int64_t current = anchors[0];
      assign_checked(q, 0, k, current);
      for (int t = 1; t < n_time; ++t) {
        current += residuals[t - 1];
        assign_checked(q, t, k, current);
      }
    } else {
      assign_checked(q, 0, k, anchors[0]);
      assign_checked(q, 1, k, anchors[1]);
      int64_t previous2 = anchors[0];
      int64_t previous1 = anchors[1];
      for (int t = 2; t < n_time; ++t) {
        const int64_t current = residuals[t - 2] + 2LL * previous1 - previous2;
        assign_checked(q, t, k, current);
        previous2 = previous1;
        previous1 = current;
      }
    }
  }
  if (reader.position() != reader.size()) stop("surplus bytes after QPP tile");
  return q;
}
