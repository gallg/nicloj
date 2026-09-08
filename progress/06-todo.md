# To do

Ordered by how much everyday work each unlocks.

## Worth doing next

**1. Keep voxels in their on-disk type.**
Everything is a `double` today, which costs 4× for `int16` and 8× for `uint8`.
The fix is to let `NdArray` hold any primitive array and convert on access, or
to add a parallel `RawArray` that `raw-data` returns and `fdata` widens. That
would also close the `int64`/`uint64` precision gap noted in
[05-status.md](05-status.md), since raw values would never pass through a
double. This is the single biggest practical limitation — it decides whether a 4D fMRI series
fits in memory. Touches `nicloj/array` and the `fdata`/`raw-data` contract, so
worth doing before anything builds on the current representation.

**2. Partial reads.**
`img.dataobj[..., 5]` in nibabel reads one volume off disk. The equivalent here
would be a `read-volume`/`read-slab` that seeks to the right byte range for
uncompressed files and streams past it for gzipped ones. Very useful for 4D
data, and cheap given the reader already computes offsets.

**3. Resampling (`nibabel.processing`).**
`resample-from-to` and `resample-to-output` need trilinear (and ideally
nearest-neighbour and spline) interpolation over the composed affine. The
geometry is already all here — `mmul`, `minverse`, `apply-affine`,
`voxel->world` — so this is mostly the interpolation kernel plus boundary
handling. `smooth-image` needs a separable Gaussian on top.

**4. Complex and colour voxel types.**
`complex64`/`complex128` show up in raw MRI work; `rgb24`/`rgba32` in
visualisation output. Both need `NdArray` to carry more than one component per
voxel, so they fall out naturally from item 1 if that redesign allows a
trailing component axis.

## Smaller, self-contained

**5. Typed header extensions.**
Give the common `ecode`s (2 DICOM, 4 AFNI, 6 comment, 18 Eprime) a parsed
representation while keeping the opaque fallback.

**6. Match nibabel's scaler derivation exactly.**
Port the `shared_range` and `floor_exact` logic from `nibabel.arraywriters` so
`:auto` picks bit-identical `scl_slope`/`scl_inter`. Only matters for
byte-comparing nicloj's output against nibabel's; values already agree.

**7. Analyze 7.5 reading.**
A 348-byte header with no NIfTI magic. The struct is already modelled — the
Analyze fields are kept for round-tripping — so this is mainly relaxing the
magic check and defaulting to the base affine.

**8. `cal_min`/`cal_max` on write.**
nibabel leaves them alone; some viewers want them set. An opt-in
`:set-cal-range? true` on `save` would compute them from the data.

**9. A `slicer` that drops axes.**
`slice-image` keeps every axis. nibabel's `img.slicer[..., 3]` drops the indexed
one. Adding integer specs alongside the `[start stop step]` triples would close
the gap.

## Infrastructure

**10. Benchmarks.**
No performance numbers exist yet. Worth measuring decode throughput against
nibabel on the MNI template before optimising anything, so item 1 can be
justified with figures rather than reasoning.

**11. Publish.**
Either a `pom.xml`/`build.clj` `jar` task for Clojars, or just a
git-dependency coordinate in the README. `deps.edn` has no licence or version
metadata yet; a `jar` task would need both.

**12. Corpus without nilearn.**
Two fixtures come from nilearn's bundled data. If nilearn moves them the
generator skips those files silently. Either vendor small crops into the repo or
make the generator fail loudly when they are missing.
