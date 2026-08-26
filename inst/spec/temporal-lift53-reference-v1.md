# Partitioned integer temporal lifting reference, version 1

Status: experimental reference for QPP-11. This document freezes the arithmetic
needed to evaluate `multiresolution-run`; it is not yet an LNA wire profile.

## Domain and partition contract

The input to one transform is a signed-int32 matrix `q[T,K]` whose rows are TRs
and columns are parcels. A file contains an ordered list of declared
independence units (normally runs). Forward transform, band quantization,
coefficient ordering, entropy coding, and inverse transform are applied to each
unit independently. No predictor, endpoint extension, lifting update, band, or
rate decision may read across a unit boundary.

Full-coefficient reference arithmetic is exact in the binary64 safe-integer
domain. A production decoder must use checked signed-int64 arithmetic and reject
overflow. Coefficients are integers; there is no floating lifting state.

The transform's lossy support is the complete declared independence unit. A
valid cross-validation fold must therefore be a union of whole transformed
units. This differs from QPP, whose only lossy operation is pointwise.

## CDF/LeGall 5/3 lifting

For one vector, split zero-based even samples into `s[i]` and odd samples into
`d[i]`. At the right boundary, duplicate the final even sample. The reversible
predict step is

```
d[i] = odd[i] - floor((s[i] + s[min(i+1,last_s)]) / 2).
```

For the update, duplicate the nearest detail at both boundaries:

```
s[i] = s[i] + floor((d[max(i-1,0)] + d[min(i,last_d)] + 2) / 4).
```

Inverse lifting subtracts the same update before adding the same prediction.
This is the reversible JPEG 2000 5/3 rounding rule. `floor` always means toward
negative infinity; language-default integer division is not an adequate
specification for negative numerators.

One level stores `coarse || detail`. Further levels recursively transform only
the coarse prefix while its length is at least two. For `L` applied levels the
final ordering is:

```
coarse-L, detail-L, detail-(L-1), ..., detail-1
```

Odd lengths retain the unmatched final even sample in the coarse band. Length
zero is forbidden. Length one is the identity. Requested levels beyond the
available dyadic depth are deterministically capped.

## Integer Haar comparator

For paired samples, the comparator uses:

```
detail = odd - even
coarse = even + floor(detail / 2)
```

The inverse subtracts the update and then adds the detail. An unmatched final
even sample passes unchanged into the coarse band. Haar uses the same recursive
ordering and partition contract as 5/3.

## Lossy and progressive semantics

The coarse band is retained with integer step one. Each detail level declares a
positive integer step. For coefficient `c` and level step `step[l]`:

```
q_c = round_ties_to_even(c / step[l])
q_c = 0 when abs(c) <= dead_zone * step[l]
c_hat = q_c * step[l]
```

The profile may declare the finest `F` detail levels dropped, in which case
their reconstructed coefficients are exactly zero. Metadata must store every
step, the dead-zone multiplier, the explicit dropped-level set, applied level
count, boundary rule, coefficient ordering, integer domain, and independence
unit ID.

Progressive decoding begins with `coarse-L` and adds complete named detail
bands. A decoder never infers a missing band from payload length. Asking for a
progressive resolution is equivalent to setting every unrequested detail band
to zero before the exact inverse. The decoded series always has `T` rows; this
is rate/resolution reduction, not an ambiguous change in logical time shape.

## Reference and required tests

`inst/prototypes/temporal-lift53-reference.R` is the executable oracle. Tests
cover lengths 1 through 17, odd and even units, zero through excessive requested
levels, negative values, multiple parcels, exact inverse, golden boundary
coefficients, coefficient coverage/order, partition isolation, band dropping,
progressive reconstruction, invalid steps, and source-domain failures.

Before wire integration, QPP-12 must add independent production code, golden
bytes, checked-int64 overflow cases, descriptor/schema validation, physical
range reads by band/unit, corruption fixtures, and installed-package parity.
