# Testing

nicloj is checked against nibabel in both directions. That is the whole point
of the harness: an independent implementation of a binary format is only as
good as its agreement with the reference.

## The two directions

**Reading.** `scripts/gen_testdata.py` uses nibabel to build a corpus in
`test-data/` and writes `test-data/manifest.edn`, recording for every file what
nibabel reports: shape, datatype, zooms, affine, qform and sform with their
codes, axis codes, `io_orientation`, units, `dim_info`, `descrip`, intent and
slice codes, extension count, the on-disk `vox_offset` and scalers, voxel
statistics, 24 sampled voxel values, and the same again after
`as_closest_canonical`. The Clojure suite loads each file and asserts nicloj
matches every one of those.

**Writing.** `nicloj.write-test` rewrites every corpus image in four layouts
and leaves them in `test-data/out/` alongside `written.edn`.
`scripts/verify_roundtrip.py` then has nibabel load each output and its
original and compare shape, datatype, affine, zooms, 25 header fields,
extension count, NIfTI version, single-file-versus-pair magic, and every voxel
at `atol=0` — exact equality.

Because the round-trip variants use `:scaling :keep`, "exact" really is exact,
not "close enough".

## Running it

```sh
./dev/check.sh
```

That runs, in order: `clojure -T:build javac`, the corpus generator if
`test-data/manifest.edn` is missing, `clojure -M:test`, then
`scripts/verify_roundtrip.py`. Override the tools with `CLOJURE=` and
`PYTHON=` if they are not on `PATH` under those names.

Regenerate the corpus deliberately after changing the generator:

```sh
python scripts/gen_testdata.py
```

The Java modules are **not** recompiled by `clojure -M:test`. After editing any
`.java` file, run `clojure -T:build javac` or use `dev/check.sh`. Forgetting
this produced two confusing test failures during development.

## The corpus

29 images, all built or copied by the generator. The synthetic volumes are
seeded from sha256 of their shape and datatype, so `gen_testdata.py` is
byte-for-byte reproducible: generating twice gives identical files and an
identical manifest.

The fixtures are:

- NIfTI-1 single file: 3D float32, gzipped int16 with `scl_slope 0.25` and
  `scl_inter -3.5`, 4D float64, 5D float32
- `.hdr`/`.img` pairs, plain and gzipped
- a big-endian NIfTI-1
- NIfTI-2 single file, 3D and 4D
- one file per supported datatype: uint8, int8, int16, uint16, int32, uint32,
  int64, uint64, float32, float64
- an oblique affine (25° about z, −12° about y, anisotropic zooms)
- qform-only, sform-only, and neither-form files
- three non-RAS orientations: LIA, PSR, ALS
- a file carrying a header extension
- two real images shipped with nilearn: the MNI152 T1 template
  (197×233×189 uint8) and `image_10426` (53×63×46 float32)

The real images matter — they are the only fixtures nicloj did not indirectly
help create, and the MNI template is large enough to catch anything that only
breaks at scale.

## The Clojure suite

`clojure -M:test` runs `test/test_runner.clj`, which is a plain
`clojure.test` runner so the project needs no test-framework dependency.

| Namespace | Covers |
| --- | --- |
| `nicloj.array-test` | Column-major layout, reshape, transpose, flip, slice, concat, the codec in both byte orders for every datatype, rounding and clamping, NaN handling, auto-scaling, error cases |
| `nicloj.affine-test` | Inverse, determinant, SVD reconstruction, polar decomposition, quaternion round trips, qform decomposition including left-handed and sheared affines, `shape_zoom_affine`, all the orientation functions |
| `nicloj.header-test` | Defaults against the NIfTI reference, purity of the accessors, qform/sform behaviour, scaling and unit fields, the binary round trip in both versions and byte orders, string truncation, extensions, rejection of malformed input and of dimensions too large for NIfTI-1 |
| `nicloj.read-test` | The manifest comparison, laziness, pair resolution from either half, rejection of unreadable input |
| `nicloj.write-test` | The four-layout round trip, `:keep` versus `:auto` versus explicit scaling, `vox_offset` and magic per layout, gzip by file name, `->bytes`, big-endian output, writing a freshly built image |
| `nicloj.ops-test` | Image construction from arrays and nested vectors, `with-*`, coordinate mapping, squeeze, concat, four-to-three, slicing with affine adjustment, reorientation preserving world positions, `dim_info` remapping |

`test/nicloj/fixtures.clj` holds the manifest reader and the approximate
comparison helpers.

## Tolerances, and why they vary

- **Voxel values from disk**: `1e-12`, essentially exact. Both sides widen the
  same stored bits to double.
- **Voxel statistics**: `1e-9` relative, with Kahan summation on the nicloj
  side. numpy sums pairwise; a naive sequential sum over the MNI template's 8.7
  million voxels drifts enough to matter, so the test compensates rather than
  loosening the bound.
- **Affines and zooms**: `1e-6` relative. NIfTI-1 stores these as float32, and
  the qform path goes through a square root and a polar decomposition.
- **Round-tripped voxels**: `0.0`. Exact, by construction, thanks to
  `:scaling :keep`.
- **Auto-scaled conversions**: one quantisation step of the target type.

## Verifying the verifier

A verifier that passes vacuously is worse than none. Two checks were run
against `scripts/verify_roundtrip.py`:

1. Its first version silently matched nothing, because the regex that pulls
   entries out of `written.edn` did not allow the commas Clojure's `pr-str`
   emits between map entries. It reported "checked 0 files"; the count in the
   output exists so that cannot pass unnoticed, and the script now also fails
   if any entry is unparsable.
2. With that fixed, one output file was corrupted by hand — one voxel and one
   `srow_x` element — and the verifier reported all three consequences (affine,
   voxels, header field) before the file was restored.
