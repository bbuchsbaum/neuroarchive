#include <Rcpp.h>
#include <cstring>
using namespace Rcpp;

// [[Rcpp::export]]
std::string morton_indices_to_hash_rcpp(IntegerVector voxelIdx) {
    // Convert integer vector to raw bytes (same layout as C int array)
    int n = voxelIdx.size();
    RawVector raw_bytes(n * sizeof(int));
    std::memcpy(&raw_bytes[0], &voxelIdx[0], n * sizeof(int));

    // Call digest::digest(raw_bytes, algo = "sha1", serialize = FALSE)
    Environment digest_ns = Environment::namespace_env("digest");
    Function digest_fn = digest_ns["digest"];
    CharacterVector result = digest_fn(raw_bytes,
                                       Named("algo") = "sha1",
                                       Named("serialize") = false);

    return "sha1:" + std::string(result[0]);
}
