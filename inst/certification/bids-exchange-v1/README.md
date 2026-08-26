# Exact installed-candidate BIDS exchange certification

The repository-root
`check-artifacts-bids-exchange-v1/certification.json` is generated only after
one source tarball is built, hashed, installed into an empty library, and used
for the complete installed court in
`inst/consumer/bids-exchange-installed-court.R`. The evidence directory is
excluded by `.Rbuildignore`, so writing the post-build receipt does not change
the candidate tarball it certifies.

The court creates its own raw and derivative fMRI sources, packs them, queries
the derivative record, performs an append update, deletes all source data,
transfers and authenticates the archive, restores a closed BIDS tree, and
checks lossless values and event bytes. The restored tree is then checked with
the official BIDS Validator 1.15.0. The final receipt embeds:

- the exact tarball SHA-256 and installed package/dependency versions;
- manifest, declaration, plan, and transfer inventory identities;
- the update and restoration receipts;
- the complete machine-readable validator result;
- zero-error and warning-key gates with a disposition for every warning; and
- the source and installed-consumer SHA-256 identities used by the court.

The certification is invalid if any step is rerun from a different source
tree or installed library. Temporary paths are working material only; the
repository-root certification receipt is the durable evidence artifact.
