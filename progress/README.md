# nicloj development notes

| File | What it covers |
| --- | --- |
| [01-architecture.md](01-architecture.md) | Module layout, why the boundaries fall where they do, build setup |
| [02-nifti-format.md](02-nifti-format.md) | The on-disk format as nicloj implements it, and every layout decision |
| [03-nibabel-parity.md](03-nibabel-parity.md) | What nibabel functionality exists here, and where behaviour deliberately differs |
| [04-testing.md](04-testing.md) | The cross-validation harness against nibabel and how to run it |
| [05-status.md](05-status.md) | What works today, verified numbers, known limitations |
| [06-todo.md](06-todo.md) | Remaining work, ordered by usefulness |
| [07-log.md](07-log.md) | Chronological record of how it was built, including the bugs found |

## The short version

nicloj reads and writes NIfTI-1 and NIfTI-2 from scratch — no external
dependency beyond Clojure itself. Byte-level work (header structs, voxel codec,
matrix and quaternion maths) is Java for speed and clarity; the API, file
handling and image operations are Clojure.

It is validated against nibabel 5.4.2 in two directions: nicloj reads a
nibabel-generated corpus and must match it field by field and voxel by voxel,
then nibabel re-reads everything nicloj writes and must find it identical to
the original. Both directions pass with exact voxel equality.
