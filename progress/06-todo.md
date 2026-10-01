# To do

Ordered by how much everyday work each unlocks.

Done: **keep voxels in their on-disk type** — `NdArray` now delegates to a
`nicloj.array.Store` per NIfTI type. See the invariant in `CLAUDE.md` and the
memory note in [05-status.md](05-status.md).

Done: **partial reads** along the last axis -- `load-slab`. Other axes still
read the whole block and slice in memory.

## Worth doing next

**1. ~~Partial reads.~~** Done for the last axis, see above.
`img.dataobj[..., 5]` in nibabel reads one volume off disk. The equivalent here
would be a `read-volume`/`read-slab` that seeks to the right byte range for
uncompressed files and streams past it for gzipped ones. Very useful for 4D
data, and cheap given the reader already computes offsets.

**2. ~~Resampling.~~** Done at `order` 0 and 1 (`nicloj/array/Resample.java`).
Cubic splines, `smooth-image` and `conform` are deliberately not planned:
linear resampling covers what nicloj is for.

**3. Complex and colour voxel types.**
`complex64`/`complex128` show up in raw MRI work; `rgb24`/`rgba32` in
visualisation output. Both need `NdArray` to carry more than one component per
voxel, so they need `Store` to carry more than one component per
voxel, which the type-per-store layout now makes straightforward.

## Smaller, self-contained

**4. Typed header extensions.**
Give the common `ecode`s (2 DICOM, 4 AFNI, 6 comment, 18 Eprime) a parsed
representation while keeping the opaque fallback.

**5. Match nibabel's scaler derivation exactly.**
Port the `shared_range` and `floor_exact` logic from `nibabel.arraywriters` so
`:auto` picks bit-identical `scl_slope`/`scl_inter`. Only matters for
byte-comparing nicloj's output against nibabel's; values already agree.

**6. ~~Analyze 7.5 reading.~~** Not planned.

**7. `cal_min`/`cal_max` on write.**
nibabel leaves them alone; some viewers want them set. An opt-in
`:set-cal-range? true` on `save` would compute them from the data.

**8. ~~A `slicer` that drops axes.~~** Done: `slice-image` takes an integer
past the third axis, as nibabel's `img.slicer[..., 3]` does.

## Infrastructure

**9. Benchmarks.**
No performance numbers exist yet. Worth measuring decode throughput against
nibabel on the MNI template before optimising anything. The one figure so far:
the 2 mm template holds 903k `int16` voxels in 1.6 MB, against 6.9 MB when
every voxel was a double.

**10. ~~Publish.~~** Not planned: nicloj stays an open repo. The version
lives in `build.clj` (`version`, now `0.1.0`), to be bumped with each release.

**11. ~~Corpus without nilearn.~~** Done: the generator now refuses to run,
before deleting anything, when nilearn's two bundled images are missing.
Skipping them had left a corpus that still passed with 373 fewer assertions.
