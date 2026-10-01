# nibabel parity

nicloj aims at the parts of nibabel used in everyday work with NIfTI volumes.
This is the mapping, and the places behaviour deliberately differs.

## What maps onto what

### Loading and saving

| nibabel | nicloj |
| --- | --- |
| `nib.load(path)` | `nii/load` |
| `nib.save(img, path)` | `nii/save` |
| `Nifti1Image(data, affine, header)` | `nii/image` |
| `img.to_bytes()` / `from_bytes()` | `nii/->bytes` / `nii/from-bytes` |
| `Nifti1Header.from_fileobj` | `nii/read-header` |
| `img.slicer[..., a:b]` (last axis) | `nii/load-slab path a b` |

### Image

| nibabel | nicloj |
| --- | --- |
| `img.shape`, `img.ndim` | `nii/shape`, `nii/ndim` |
| `img.affine` | `nii/affine` |
| `img.header` | `nii/header` |
| `img.get_fdata()` | `nii/fdata` |
| `np.asarray(img.dataobj)`, `dataobj.get_unscaled()` | `nii/raw-data` |
| `img.get_data_dtype()` / `set_data_dtype` | `nii/data-dtype` / `nii/with-dtype` |
| `img.header.get_zooms()` | `nii/zooms` |
| — | `nii/voxel`, `nii/describe`, `nii/loaded?` |

### Header

| nibabel | nicloj |
| --- | --- |
| `get_qform` / `set_qform` | `nii/qform` / `nii/set-qform` |
| `get_sform` / `set_sform` | `nii/sform` / `nii/set-sform` |
| `get_best_affine` | `nii/best-affine` |
| `get_base_affine` | `nii/base-affine` |
| `get_zooms` / `set_zooms` | `nii/set-zooms` |
| `get_data_shape` / `set_data_shape` | `nii/set-shape` |
| `get_slope_inter` / `set_slope_inter` | `nii/slope-inter` / `nii/set-slope-inter` |
| `get_xyzt_units` / `set_xyzt_units` | `nii/xyzt-units` / `nii/set-xyzt-units` |
| `get_dim_info` / `set_dim_info` | `nii/dim-info` / `nii/set-dim-info` |
| `hdr.extensions` | `nii/extensions` |
| `print(hdr)` | `nii/describe-header`, `nii/header->map` |

### Orientation (`nibabel.orientations`)

| nibabel | nicloj |
| --- | --- |
| `io_orientation` | `nii/io-orientation` |
| `aff2axcodes` | `nii/axcodes` |
| `ornt2axcodes` / `axcodes2ornt` | `nii/ornt->axcodes` / `nii/axcodes->ornt` |
| `ornt_transform` | `nii/ornt-transform` |
| `inv_ornt_aff` | `nii/inv-ornt-aff` |
| `apply_orientation` | `nii/apply-orientation` |
| `img.as_reoriented` | `nii/as-reoriented` |
| `nib.as_closest_canonical` | `nii/as-closest-canonical` |

### Whole-image operations (`nibabel.funcs`)

| nibabel | nicloj |
| --- | --- |
| `squeeze_image` | `nii/squeeze-image` |
| `concat_images` | `nii/concat-images` |
| `four_to_three` | `nii/four-to-three` |
| `img.slicer[...]`, `img.slicer[..., 3]` | `nii/slice-image`, with an integer to drop an axis past the third |

### Resampling (`nibabel.processing`)

| nibabel | nicloj |
| --- | --- |
| `resample_from_to(img, to, order, cval=...)` | `nii/resample-from-to img to :order :cval` |
| `resample_to_output(img, voxel_sizes, order)` | `nii/resample-to-output img :voxel-sizes :order` |

Only `order` 0 and 1 exist, and the default is 1 where nibabel's is 3. At
those orders the voxels match nibabel's to rounding noise: exactly for integer
types, within 1e-12 for floats, over the corpus and 300 random cases.

### Affines (`nibabel.affines`)

| nibabel | nicloj |
| --- | --- |
| `apply_affine` | `nii/apply-affine` |
| `from_matvec` / `to_matvec` | `nii/from-mat-vec` / `nii/to-mat-vec` |
| `voxel_sizes` | `nii/voxel-sizes` |
| `npl.inv(aff)`, `aff @ bff` | `nii/minverse`, `nii/mmul` |
| — | `nii/voxel->world`, `nii/world->voxel` |

## Deliberate differences

**Scaling on load.** nibabel moves `scl_slope`/`scl_inter` into its array proxy
and blanks the header fields it hands you, so `hdr.get_slope_inter()` returns
`None` even for a file that clearly has scaling on disk. nicloj leaves the
fields as they were read, so `nii/slope-inter` reports what the file says.
`fdata` applies them either way, so voxel values agree exactly. This tripped up
the first version of the test suite, which compared against nibabel's blanked
fields; the corpus manifest now records the on-disk bytes instead.

**`vox_offset` on load.** Same story: nibabel resets it to 0 after reading a
single-file image. nicloj reports the real offset.

**Scaling on save.** nibabel derives `scl_slope`/`scl_inter` through its
`ArrayWriter` machinery, which works hard to make the scalers exactly
representable in the header's float32. nicloj's `:auto` computes the same
linear map but relies on clamping in the encoder instead of chasing exact
representability, so a float-to-integer conversion can differ from nibabel's in
the last bit or two of the scalers. Values still round-trip within one
quantisation step. Where an exact rewrite matters, `:scaling :keep` writes the
stored voxels back untouched and reuses the header's scalers — the round-trip
tests use it and get byte-identical voxels.

**Explicit scaling control.** nicloj's `save` takes `:scaling`, which nibabel
has no direct equivalent for: `:auto` (default), `:keep`, or an explicit
`[slope inter]` pair.

**`scl_slope` when no scaling is wanted.** nibabel writes `NaN`; nicloj writes
`1.0` with `scl_inter` `0.0`. Both mean "no scaling" under the spec, and both
read back identically, but the bytes differ.

**Value type.** Voxels are stored in the file's own type, as in nibabel, but
every accessor hands back a `double`. So `nd/dtype` answers what
`get_data_dtype()` does, while there is no equivalent of asking for the values
themselves in another width, the way `get_fdata(dtype=np.float32)` does.

**xform codes after slicing.** nibabel's slicer rebuilds the image from its
affine, so a qform-only file comes back with `qform_code` 0 and an `:aligned`
sform. `load-slab` shifts whichever forms were set and keeps their codes. The
affines agree; only the codes differ.

**`dim[0] = 0`.** nibabel reads it as shape `(0,)`, zero voxels; nicloj reads
it as a 0-d image of one voxel, which is also what a fresh `new-header`
describes. Changing it would change every new header for an edge case where
nibabel itself drops the voxel nicloj keeps.

**Resampling at the edge.** scipy gives `cval` to a point even 1e-15 voxels
past the grid, so when composing affines leaves rounding residue it can lose a
whole edge slice; resampling an image onto its own grid did exactly that.
nicloj counts a point within 1e-9 voxels of the edge as on it.

**Resampling a series onto a 3D grid.** nibabel needs the target to have as
many axes as the source; nicloj fills missing trailing axes in from the
source, so a 4D series resamples onto a 3D template volume by volume.

**Freesurfer hacks.** nibabel special-cases two Freesurfer conventions in
`get_data_shape`/`set_data_shape`: `dim[1] == -1` with the real length in
`glmin`, and the ico7 surface shape `(27307, 1, 6)` standing for
`(163842, 1, 1)`. nicloj does not implement either; such files read with their
literal `dim`.

**No image classes beyond NIfTI.** nibabel also reads MINC, ECAT, PAR/REC,
GIFTI, CIFTI-2, Analyze and Freesurfer formats. nicloj is NIfTI only.

## Not implemented

- Complex (`complex64`, `complex128`) and colour (`rgb24`, `rgba32`) voxel
  types — recognised, but decoding raises.
- `nibabel.processing`: cubic-spline resampling (`order` 2 to 5),
  `smooth_image` and `conform`.
- Typed extension classes (`Nifti1Extension` subclasses for AFNI, DICOM,
  comments). Extensions are carried as opaque bytes with their `ecode`.
- `nibabel.imagestats`, `nibabel.viewers`, the command-line tools.
- Memory-mapped and partially-read (`fileslice`) access. `load` is lazy about
  *when* it reads the voxel block, but reads all of it when it does.
