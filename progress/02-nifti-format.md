# The NIfTI format as nicloj implements it

Field names, offsets and semantics come from the public-domain `nifti1.h` and
`nifti2.h` reference headers at <https://nifti.nimh.nih.gov/>. This file records
what nicloj does with them, so the next person does not have to re-derive it.

## Detecting the layout

`sizeof_hdr`, the first four bytes, says everything:

| Value read little-endian | Meaning |
| --- | --- |
| `348` | NIfTI-1, little-endian |
| `540` | NIfTI-2, little-endian |
| `0x5C010000` (348 byte-swapped) | NIfTI-1, big-endian |
| `0x1C020000` (540 byte-swapped) | NIfTI-2, big-endian |

Anything else is not a NIfTI file. `NiftiHeader.probe` does this; the reader
uses it to decide how many more bytes to pull before parsing.

## Field layout

`NiftiHeader` holds the union of both versions using the wider representation
(`long`, `double`) and narrows on write. NIfTI-1's Analyze 7.5 leftovers
(`data_type`, `db_name`, `extents`, `session_error`, `regular`, `glmax`,
`glmin`) are kept as fields so they survive a round trip rather than being
silently zeroed.

**NIfTI-1, 348 bytes.** `sizeof_hdr`(0) `data_type`(4) `db_name`(14)
`extents`(32) `session_error`(36) `regular`(38) `dim_info`(39) `dim[8]`(40, int16)
`intent_p1..3`(56, float32) `intent_code`(68) `datatype`(70) `bitpix`(72)
`slice_start`(74) `pixdim[8]`(76, float32) `vox_offset`(108, float32)
`scl_slope`(112) `scl_inter`(116) `slice_end`(120) `slice_code`(122)
`xyzt_units`(123) `cal_max`(124) `cal_min`(128) `slice_duration`(132)
`toffset`(136) `glmax`(140) `glmin`(144) `descrip[80]`(148) `aux_file[24]`(228)
`qform_code`(252) `sform_code`(254) `quatern_b/c/d`(256) `qoffset_x/y/z`(268)
`srow_x/y/z[4]`(280) `intent_name[16]`(328) `magic[4]`(344).

**NIfTI-2, 540 bytes.** Same information, reordered and widened:
`sizeof_hdr`(0) `magic[8]`(4) `datatype`(12) `bitpix`(14) `dim[8]`(16, int64)
`intent_p1..3`(80, float64) `pixdim[8]`(104, float64) `vox_offset`(168, int64)
`scl_slope`(176) `scl_inter`(184) `cal_max`(192) `cal_min`(200)
`slice_duration`(208) `toffset`(216) `slice_start`(224) `slice_end`(232)
`descrip[80]`(240) `aux_file[24]`(320) `qform_code`(344) `sform_code`(348)
`quatern_b/c/d`(352) `qoffset_x/y/z`(376) `srow_x/y/z[4]`(400)
`slice_code`(496) `xyzt_units`(500) `intent_code`(504) `intent_name[16]`(508)
`dim_info`(524) `unused_str[15]`(525).

Both were checked against nibabel's own structured dtypes rather than trusted
from memory.

Notable per-version quirks nicloj has to respect:

- NIfTI-1 `dim` is `int16`, so a dimension over 32767 cannot be written; the
  writer raises a clear error naming NIfTI-2 as the fix.
- NIfTI-1 `vox_offset` is a **float32**, which is why offsets are exact only up
  to 2^24. NIfTI-2 makes it an `int64`.
- `dim_info`, `slice_code` and `xyzt_units` are unsigned bytes in NIfTI-1 and
  `int32` in NIfTI-2; they are read masked with `0xFF` in the former.
- NIfTI-2's magic is 8 bytes: the 4-character tag plus `0D 0A 1A 0A`, a
  CR/LF/EOF guard that makes a text-mode transfer corrupt the file detectably.

String fields are fixed-length, NUL-terminated, and decoded as Latin-1 so
arbitrary bytes survive; writing truncates rather than overflowing.

## File layouts

| Name | `magic` | Where voxels live |
| --- | --- | --- |
| `.nii`, `.nii.gz` | `n+1` / `n+2` | same file, at `vox_offset` |
| `.hdr` + `.img` | `ni1` / `ni2` | the `.img` file, at `vox_offset` (normally 0) |

For a single file, `vox_offset` is at least `sizeof_hdr + 4`: 352 for NIfTI-1,
544 for NIfTI-2, more when extensions are present. The reader tolerates a
`vox_offset` below that minimum by starting the data right after the 4-byte
extender, which is where nibabel's fix-up puts it too (nibabel then refuses
the file; nicloj reads it). A negative `vox_offset` in a pair is an error.

Trailing NUL bytes of an extension's content are padding, and are stripped on
read as nibabel does; the writer pads again to the 16-byte boundary.

Compression is detected by sniffing the gzip magic `1F 8B` rather than trusting
the file name, so a mislabelled `.nii` holding compressed bytes still loads.
Writing compresses when the name ends in `.gz`. Because gzip streams are not
seekable, reading is always sequential: skip to the offset, then read.

## Extensions

After the fixed header comes a four-byte extender. A non-zero first byte means
extension blocks follow, each `esize:int32, ecode:int32, payload`, padded to a
multiple of 16 bytes. nicloj parses them into `Extension` records, keeps the
payload opaque, and writes them back unchanged — including recomputing
`vox_offset` to make room. In a `.hdr`/`.img` pair, the extension region is
simply the rest of the `.hdr` file.

## Datatypes

Supported: `uint8`(2) `int16`(4) `int32`(8) `float32`(16) `float64`(64)
`int8`(256) `uint16`(512) `uint32`(768) `int64`(1024) `uint64`(1280).

Recognised but not yet decodable: `complex64`(32) `rgb24`(128)
`complex128`(1792) `rgba32`(2304) — these raise a clear error.

Not in the table at all: `float128`(1536) and `complex256`(2048), which have no
Java equivalent, and `binary`(1), which nibabel does not support for NIfTI
either.

Unsigned types are widened through masking (`& 0xFF`, `& 0xFFFF`,
`& 0xFFFFFFFFL`); `uint64` adds 2^64 to negative longs. `uint64`'s maximum is
reported as the largest double below 2^64 rather than 2^64−1, because 2^64−1 is
not representable and rounding up to 2^64 made the encoder wrap to zero. That
was a real bug caught by the datatype round-trip test.

Writing an integer type rounds half-to-even (matching numpy's `rint`, which
nibabel uses) and clamps to the type range; NaN becomes zero.

## Affines

Three sources, in the order `best-affine` prefers them:

1. **sform** — `srow_x/y/z` read straight off as the top three rows of a 4×4.
   Used whenever `sform_code != 0`.
2. **qform** — a unit quaternion, a translation and `pixdim`. Used when
   `sform_code == 0` and `qform_code != 0`. The rotation is
   `quat2mat(fillpositive(b, c, d))`, its columns scaled by
   `pixdim[1]`, `pixdim[2]` and `qfac * pixdim[3]`, where `qfac = pixdim[0]`
   carries the handedness. Only `b`, `c`, `d` are stored, so the real part is
   recovered as `sqrt(1 - (b² + c² + d²))`, with nibabel's `3.5762787e-07`
   threshold absorbing float32 storage error.
3. **base affine** — a diagonal affine centred on the image with the Analyze
   left-right flip, when both codes are 0.

Going the other way (`Affines.toQform`), an affine is decomposed into column
norms (the zooms) and a normalised rotation; if the determinant is negative,
`qfac` becomes −1 and the third column is negated. The rotation is then
projected onto the nearest orthogonal matrix, because the qform cannot express
shear — the result reports whether any shear was dropped.

Writing an image the way nicloj does it puts the affine in both places: sform
with code `:aligned`, qform values with code `:unknown`. This matches nibabel's
`_affine2header`.

## Linear algebra written for this

- **Thin SVD** by one-sided Jacobi rotations (`Svd`). Needed for the polar
  decomposition in both `toQform` and `io_orientation`, including the
  rank-deficient case that a Newton polar iteration cannot handle. Transposes
  the input when it has more columns than rows.
- **Quaternion from matrix** by Shepperd's method, with the real part
  normalised non-negative. nibabel uses an eigendecomposition of a symmetric
  4×4; for a proper rotation matrix — which is all either code ever sees, since
  both orthogonalise first — the two agree exactly.
- **Inverse and determinant** by Gauss-Jordan and LU, both with partial
  pivoting.

## Orientations

An orientation is one `[out-axis flip]` row per voxel axis. `io_orientation`
normalises the affine's rotation block by its column norms, takes the polar
factor, and then claims output axes greedily from the strongest input axis
down, zeroing each claimed output row so ties resolve consistently. The
strongest-first ordering is not incidental: nibabel added it precisely so that
an axis gets the same label regardless of the order the axes happen to sit in,
and dropping it would make nicloj disagree on affines with near-ties.

`apply-orientation` flips first, then transposes by the inverse permutation,
and carries any axes past the third along untouched — so it works on 4D and 5D
images. `inv-ornt-aff` produces the affine that undoes those moves, so
`new-affine = old-affine × inv-ornt-aff(ornt, shape)` keeps every voxel at the
same world coordinate.
