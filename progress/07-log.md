# Development log

A record of the order things were built and what went wrong, so the reasoning
behind the current shape is recoverable.

## 1. Reference gathering

Before writing anything, the authoritative behaviour was read out of the
installed nibabel 5.4.2 rather than recalled:

- `nifti1.py` — `get_qform`, `set_qform`, `get_sform`, `set_sform`,
  `get_best_affine`, `get_slope_inter`, `get_xyzt_units`, `get_zooms`,
  `get_dim_info`, `_affine2header`, `default_structarr`
- `quaternions.py` — `fillpositive`, `quat2mat`, `mat2quat`
- `affines.py` — `from_matvec`, `to_matvec`, `apply_affine`
- `analyze.py` — `shape_zoom_affine`, `set_data_dtype`
- `orientations.py` — all of it
- `arraywriters.py` — `_do_scaling`, `_range_scale`, `_iu2iu`
- both header dtypes, dumped field by field with their offsets
- the code-label vocabularies for xform, units, slice order and intent, which
  turned out to match `Codes.java` exactly

nibabel's licence was confirmed as MIT from its package metadata (there is no
`LICENSE` file in the wheel), which is what makes reimplementing these
algorithms with attribution acceptable. `nifti1.h`/`nifti2.h` are public domain.

nilearn's bundled `datasets/data/` directory was checked for real NIfTI files
usable as fixtures without a download; it has four, two of which are now in the
corpus.

## 2. Java modules

Written bottom-up: `DataType` and `NiftiError`, then `Extension` and `Codes`,
then `NiftiHeader`; then `NdArray` and `Codec`; then `Mat`, `Svd`,
`Quaternions`, `Affines` and `Orientations`.

Two decisions taken here shaped everything after:

- `NiftiHeader` is one mutable class covering both versions, with public
  fields, rather than a class hierarchy. Clojure's `set!` works on public
  fields, so the functional wrapper in `nicloj.core.header` is thin.
- `Svd` implements one-sided Jacobi rather than a Newton polar iteration.
  Newton is shorter, but `io_orientation` needs singular values thresholded to
  detect dropped axes, and a rank-deficient matrix breaks the Newton iteration.

`Quaternions.fromMatrix` uses Shepperd's method where nibabel uses an
eigendecomposition. That is safe only because both sides orthogonalise the
matrix first, so the input is always a proper rotation — noted in the source.

## 3. Clojure modules

`io/paths` and `io/stream`, then `core/linalg`, `core/header`, `core/ndarray`,
`core/image`, `core/nifti`, `core/orientation`, `core/funcs`, then the
`api/image` facade.

Design points settled while writing these:

- **The facade re-exports through a macro.** 60-odd names with their docstrings
  and arglists carried over, in 77 lines, instead of hand-written `def`s.
- **`save` takes an explicit `:scaling`.** The alternative — inferring whether
  to preserve or recompute the scalers — was rejected as too magical to
  document honestly. `:auto`, `:keep` and an explicit pair each mean one thing.
- **`image` and `with-data` reset `scl_slope`/`scl_inter`.** Because `fdata`
  always applies header scaling, data handed in directly has to be treated as
  already scaled or it would be scaled twice. This invariant was arrived at
  while writing `as-reoriented`, which builds a new image out of `fdata`.
- **`apply-orientation` lives in Clojure, not `nicloj/affine`.** It is the only
  orientation operation touching voxels, and putting it in Java would make the
  pure-maths module depend on `NdArray`.

## 4. Test corpus

`scripts/gen_testdata.py` builds 29 images and a manifest of what nibabel
reports for each. Getting the manifest right took three iterations:

1. **Byte order in `disk_fields`.** The raw header fields were unpacked as
   little-endian, so the big-endian fixture reported
   `scl_slope 4.6e-41` — a denormal, the tell-tale of a byte-swapped float.
   Fixed by probing `sizeof_hdr` for the order first.
2. **`vox_offset` and the scalers.** nibabel blanks both after loading a
   single-file image (it moves scaling into its array proxy), so the manifest
   was recording 0 and `NaN` for files that plainly had 352 and `0.25`/`-3.5`
   on disk. The generator now reads those four fields out of the raw bytes.
3. **`get_slope_inter()` is not usable as an expectation** for the same reason.
   The `:slope`/`:inter` entries were dropped; the on-disk `:scl` pair is the
   real check, and the voxel-value assertions already prove the scaling is
   applied.

## 5. First full run

75 tests, 3799 assertions, 10 failures. All of them informative:

- **`:scaling :keep` wrote the wrong voxels.** `resolve-scaling` returned the
  raw array and the header's slope/inter, and `encode` then *inverted* the
  scaling on data that was already raw. Fixed by having `resolve-scaling`
  return a `:raw?` flag, so `encode` writes raw arrays with identity scaling
  while still recording the scalers in the header. This is exactly the bug a
  round-trip test exists to catch.
- **`uint64` did not round-trip.** `DataType.maxValue()` returned
  `18446744073709551615.0`, which rounds to 2^64 as a double; clamping to it
  and subtracting 2^64 in the encoder produced 0. Now returns
  `Math.nextDown(0x1p64)`.
- **Two stale-class failures.** `DataType.java` and `Codec.java` had been
  edited after the last `javac`, and `clojure -M:test` does not recompile.
  `dev/check.sh` now always runs `javac` first.
- **Three test expectations were wrong, not the code.** `float32` probe values
  of `±1.5e10` are not exactly representable; `[-1.0 0.0 1.0]` is integral and
  in range so `autoScale` correctly declines to scale it; and NIfTI-1 stores
  `scl_slope` as float32, so `0.01` reads back as `0.009999999776482582`.

## 6. Second run and verification

All 3774 assertions passed, and `verify_roundtrip.py` reported "checked 0 files
— all match". Zero files is not a pass. The regex pulling entries out of
`written.edn` did not allow the commas Clojure's `pr-str` puts between map
entries. After fixing it: 116 files, all matching.

To confirm the verifier was not passing vacuously, one output file was
corrupted by hand — one voxel value and one `srow_x` element — and it reported
all three consequences before the file was restored.

## 7. Polish

Reflection warnings were checked with `*warn-on-reflection*` and the four
found were fixed: a `GZIPOutputStream` constructor needing an `int` hint, a
`double[]` row from `aget`, and the two `Extension` field reads in
`header/->map`. `to-mat-vec` was added for symmetry with `from-mat-vec`, with a
round-trip test.

## 8. Bug hunt

A deliberate pass over paths the suite did not reach. Five real defects, each
fixed with the smallest change that closed it, and each now covered by a test.
Behaviour was probed before and after every fix, and the whole suite plus the
nibabel verification re-run to check for regressions.

**`from-bytes` silently truncated data.** `java.util.Arrays.copyOfRange` pads
with zeros when the end index runs past the array, so a short buffer produced a
full-shaped image whose tail voxels were 0.0 instead of an error. Reading a
2368-byte image from a 2268-byte buffer gave `100.09` as `0.0`. Now bounds-checked
before slicing. This was the worst of the five: silent corruption, not a crash.

**`Codec.encode` overflowed its allocation size.** `n * type.itemSize()` was int
arithmetic, so a float64 array of more than 268435456 voxels wrapped negative —
`268435456 * 8` is exactly 2^31. `ByteBuffer.allocate` would then throw
`IllegalArgumentException: capacity < 0`, or for other sizes allocate too small
and fail later with `BufferOverflowException`. Now computed as a `long` and
rejected with the same message the reader already used. Verified by arithmetic
rather than by allocating 2 GB, so it has no direct test; the analogous check on
the decode side is tested.

**`existing-image` only resolved half of a mismatched pair.** Its docstring
promised tolerance of a `.hdr`/`.img` pair that disagrees about gzip, but the
fallback was applied to the image half only. A plain `x.hdr` beside a gzipped
`x.img.gz`, loaded by the image name, resolved the header to a non-existent
`x.hdr.gz` and died with `FileNotFoundException`. Replaced the one-sided logic
with an `on-disk` helper applied to both halves — which also made the function
shorter. Single-file names are deliberately left alone, so `load "x.nii"` never
silently picks up `x.nii.gz`.

**`with-header` silently rescaled voxels.** Voxel data is stored raw with
`scl_slope`/`scl_inter` in the header, so replacing the header with one carrying
different scaling changed what `fdata` returned. On the scaled int16 fixture,
voxel `[0 0 0]` went from `-3.25` to `1.0`. The old scalers now carry over,
because the stored data is relative to them; every other field of the new header
is still adopted. Laziness and the purity of the header argument were both
re-checked afterwards.

**`base-affine` crashed on a dimensionless header.** `zooms()` returns `(1.0)`
when `dim[0]` is 0 — matching nibabel — while `shape()` returns nothing, so
`shapeZoomAffine` rejected the mismatched lengths and `(nii/base-affine
(nii/new-header))` failed with "shape and zooms must have the same length".
Zooms are now trimmed to the shape length, which only ever affects this one
degenerate case.

Also probed and found sound: NIfTI-2 pairs, big-endian NIfTI-2, 7-dimensional
images, zero-voxel images, saving via an `.img` name, `slice-image` with fewer
specs than axes, and `dim[0]` inconsistent with the data. Two rough error
messages were left alone as caller errors that already fail loudly:
`nd/slice` with too few specs, and `concat-images` with an out-of-range axis.

## 9. Packaging, and what a clean clone found

Machine-specific files were moved out of the repo: `dev/in-dev.sh` (this
machine's distrobox wrapper) joined `.venv/`, `target/`, `.cpcache/` and
`test-data/` in `.gitignore`, leaving 47 tracked files — source, tests, scripts,
docs, `deps.edn`, `build.clj` and the portable `dev/check.sh`. CLAUDE.md now
carries the contents of `in-dev.sh` so it can be recreated. The README was cut
from 92 lines to 61.

To check that nothing tracked depended on something ignored, the non-ignored
files were copied to a scratch directory and `dev/check.sh` run there. It
failed — five assertions, all on `dtype-uint64.nii`. Two real defects behind it:

**The test corpus was not reproducible.** `gen_testdata.py` seeded its RNG from
`hash(str(shape) + str(dtype))`, and Python salts `hash()` per process. Every
run produced different voxel data, so the suite was quietly non-deterministic —
exactly the setup that yields "works on my machine". Now seeded from sha256:
generating the corpus twice, or in a fresh clone, gives byte-identical files and
an identical manifest, which was verified by checksum.

**Unsigned fixtures wrapped around.** `ramp(shape, uint64, lo=0, hi=200)` added
Gaussian noise to a ramp starting at 0, so values near zero could go negative
and the unsigned cast wrapped them to near 2^64 — not a deliberate boundary
test, just an accident that the old random seed usually hid. Values are now
clipped to the target type's range before casting.

That accident did expose a genuine limitation, now recorded in
[05-status.md](05-status.md) and pinned by a test: because every voxel is a
double, `int64`/`uint64` magnitudes past 2^53 are rounded on read and cannot be
written back bit for bit. nibabel's `get_fdata()` loses the same precision, but
its rewrite path copies the stored bytes, so this is a gap only in nicloj's
write path. It closes with todo item 1.

Fixing the corpus then broke two of my own tests, which had hard-coded numbers
copied out of the old data — a voxel value and a data range used as a
tolerance. Both were rewritten to assert relationships or derive the number
from the image, which is what they should have done in the first place.

Final state: 80 tests, 3789 assertions, and the nibabel verification of all 116
written files, passing both in the working tree and in a clean clone whose
corpus was generated independently.

## 10. Licence

MIT, in a `LICENSE` file whose text matches the canonical OSI wording exactly
so GitHub's licence detection recognises it. The choice was unconstrained:
`nifti1.h`/`nifti2.h` are public domain, and the algorithms reimplemented from
nibabel are MIT, so MIT-to-MIT adds no obligation beyond the attribution
already in the source comments.

The provenance note lives in the README rather than in `LICENSE`, because extra
prose in that file stops automated licence detectors matching it.
