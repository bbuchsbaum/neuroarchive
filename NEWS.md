## neuroarchive 0.1.0

* Implemented the "rock-solid quant" transform for integer quantisation.
* Quantized datasets are stored using 8- or 16-bit integer types
  with the bit depth recorded in a `quant_bits` attribute.
* Added per-voxel mode (`scale_scope = "voxel"`) which computes and
  stores scale and offset per voxel using block-wise processing.
* Clipping logic now counts clipped samples and can warn or abort
  based on configured thresholds.
* A detailed quantization report is produced in JSON format and
  written to the file (compressed with GZIP). Retrieve it with
  `lna_get_quant_report()`.
* Added `lna.debug.temporal` option. When `TRUE`, `forward_step.temporal` and
  `invert_step.temporal` emit debug messages; defaults to `FALSE` for quiet runs.
* `forward_step.spat.haar_octwave` now stores a list of valid finest-level
  block Morton codes at `/aux_meta/haar_octwave/valid_blocks_L-1` for
  optimized ROI streaming.
