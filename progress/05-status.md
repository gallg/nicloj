# Status

As of the initial implementation, the bug-fix pass, and the voxel-storage
pass in [07-log.md](07-log.md).

## Verified

```
clojure -M:test                     83 tests, 3833 assertions, 0 failures, 0 errors
python scripts/verify_roundtrip.py  116 files nicloj wrote, all match nibabel
```

`clojure -M -e "(set! *warn-on-reflection* true) (require 'nicloj.api.image)"`
is clean — no reflective call sites.

The read direction matches nibabel 5.4.2 on every field of all 29 corpus
images. The write direction produces files nibabel reads back with exactly
equal voxels and matching header fields, in four layouts each.

## What works

**Reading.** NIfTI-1 (348 bytes) and NIfTI-2 (540 bytes), little- and
big-endian, single `.nii` files and `.hdr`/`.img` pairs, gzipped or not, with
either half of a pair naming the image. Version and byte order are detected
from `sizeof_hdr`; compression from the gzip magic. Header extensions are
parsed and preserved. Voxel data is read on first use.

**Writing.** All of the above, chosen by file name and options, with header
extensions and `vox_offset` handled correctly. `:scaling :keep` reproduces a
loaded file's voxels byte for byte; `:auto` derives scalers that fit a narrower
datatype; an explicit `[slope inter]` pair is also accepted.

**Datatypes.** uint8, int8, int16, uint16, int32, uint32, int64, uint64,
float32, float64 — read and written in both byte orders.

**Geometry.** qform and sform, their codes, the centred Analyze fallback,
quaternion decomposition with handedness and shear detection, voxel-to-world
and world-to-voxel mapping, and the full orientation toolkit including
`as-closest-canonical` on images of any dimensionality.

**Operations.** `squeeze-image`, `concat-images`, `four-to-three`,
`slice-image` (cropping and subsampling with the affine adjusted so voxels keep
their world positions), and the array primitives underneath them.

## Size

| | Files | Lines |
| --- | --- | --- |
| Java modules (`header`, `array`, `affine`) | 12 | 1842 |
| Clojure modules (`io`, `core`, `api`) | 10 | 1115 |
| Tests | 7 | 954 |
| Python harness | 2 | ~400 |

Largest single file is `NiftiHeader.java` at 457 lines, most of which is the
two binary layouts written out field by field. `nicloj/api/image.clj` is 77
lines for the whole public surface, because it re-exports through a small macro
rather than restating docstrings.

## Known limitations

**Memory.** Voxels are held at their on-disk width, so a loaded image costs
about what the file costs: the 2 mm MNI152 template (903k `int16` voxels) is
1.6 MB rather than the 6.9 MB it took when everything was a `double`. Reading
an element still gives a `double`. Two things still widen: `fdata` on an image
with `scl_slope`/`scl_inter`, since scaled values no longer fit the source
type, and `nd/values`, which realises a boxed Clojure vector. A single array is
capped at `Integer.MAX_VALUE` elements, and the reader refuses a voxel block
over 2 GiB with a clear message rather than failing obscurely.

**64-bit integers past 2^53 read approximately, but round-trip exactly.**
`nd/value` and friends hand back a `double`, so a `uint64` of 2^64-2 reads as
2^64. The bits themselves are kept, and rewriting the image to the same type
copies them through `Store.getLong`, so the file survives unchanged. Converting
to a *different* integer type still goes through a double and is still lossy.
Pinned by `large-64-bit-values-survive-a-round-trip` in `nicloj.array-test`.

**No complex or colour voxels.** `complex64`, `complex128`, `rgb24` and
`rgba32` are recognised in the datatype table so headers describing them parse,
but decoding raises. `float128` and `complex256` have no Java equivalent and are
not in the table at all.

**No resampling.** Nothing in nicloj interpolates. `slice-image` crops and
subsamples on the existing grid; there is no equivalent of
`nibabel.processing.resample_from_to`.

**Auto-derived scalers are not bit-identical to nibabel's.** See
[03-nibabel-parity.md](03-nibabel-parity.md). Values agree to within a
quantisation step; the scalers themselves can differ in their last bits.

**Extensions are opaque.** Carried through with their `ecode`, but not parsed
into the typed forms nibabel offers for AFNI and DICOM payloads.

**No Freesurfer `dim` conventions.** Files using the `dim[1] == -1` large-vector
hack or the ico7 surface shape read with their literal `dim`.

**Reading is all-or-nothing.** Laziness controls *when* the voxel block is
read, not how much; there is no partial or memory-mapped access.
