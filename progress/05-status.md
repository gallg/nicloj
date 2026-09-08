# Status

As of the initial implementation plus the bug-fix pass in [07-log.md](07-log.md).

## Verified

```
clojure -M:test                     80 tests, 3789 assertions, 0 failures, 0 errors
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

**Memory.** Every voxel becomes a `double`, so an `int16` volume costs four
times its file size in RAM and a `uint8` one eight times. The MNI152 template
(8.7 million voxels) is about 70 MB loaded. A single array is also capped at
`Integer.MAX_VALUE` elements, and the reader refuses a voxel block over 2 GiB
with a clear message rather than failing obscurely.

**64-bit integers past 2^53 are approximate.** A consequence of the same
choice: `int64`/`uint64` values larger than a double can hold exactly are
rounded on read and cannot be written back bit for bit. nibabel's
`get_fdata()` loses the same precision, but its rewrite path copies the stored
bytes, so this is a gap only in nicloj's write path. Pinned by
`large-64-bit-values-are-approximate` in `nicloj.array-test`.

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
