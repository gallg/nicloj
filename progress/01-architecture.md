# Architecture

## Modules

A module is a folder under `nicloj/` holding one self-contained feature, written
entirely in Java or entirely in Clojure.

| Module | Language | Namespace / package | Responsibility |
| --- | --- | --- | --- |
| `nicloj/header` | Java | `nicloj.header` | The NIfTI-1 and NIfTI-2 header structs, the datatype table, header extensions, and the symbolic names for enumerated fields |
| `nicloj/array` | Java | `nicloj.array` | A dense n-dimensional voxel array, and conversion between raw bytes and that array |
| `nicloj/affine` | Java | `nicloj.affine` | Matrix algebra, thin SVD, quaternions, qform/sform conversion, axis orientations |
| `nicloj/io` | Clojure | `nicloj.io` | Which files make up an image, and gzip-aware byte streams |
| `nicloj/core` | Clojure | `nicloj.core` | The header API, the image type, the reader and writer, orientation and whole-image operations |
| `nicloj/api` | Clojure | `nicloj.api` | One namespace re-exporting the whole public surface |

Dependencies run strictly downward: `api` → `core` → `io`, `affine`, `array` →
`header`. The three Java modules know nothing about Clojure and can be used
from plain Java.

### Why the split falls here

**Java for the byte and number work.** A NIfTI header is a fixed C struct;
expressing it as a mutable class with public fields and a `ByteBuffer` codec is
both shorter and faster than any Clojure equivalent. The same goes for decoding
millions of voxels and for 3×3 matrix arithmetic, where boxing would dominate.

**Clojure for everything with a decision in it.** File-name resolution, the
functional header API, laziness, and image operations all read better as
Clojure, and that is where a user actually works.

**`nicloj/array` is separate from `nicloj/header`.** The codec needs the
datatype table, but the array type is useful on its own and the header struct
should not depend on voxel storage.

**`nicloj/affine` is pure maths.** It deliberately does not reference
`NdArray`, so `apply-orientation` — the one orientation operation that touches
voxels — lives in `nicloj.core.orientation` and composes `NdArray/flip` and
`NdArray/transpose` instead.

## Classpath and build

The project root is the classpath root, so the `nicloj/` source folder doubles
as the top-level package:

```
nicloj/core/nifti.clj             -> namespace nicloj.core.nifti
nicloj/header/NiftiHeader.java    -> class nicloj.header.NiftiHeader
```

`deps.edn` therefore has `:paths ["." "target/classes"]`. Java sources are
compiled ahead of time into `target/classes`, which the Clojure code then loads
classes from:

```sh
clojure -T:build javac
```

`build.clj` uses `tools.build`, but the task is a single `javac` invocation and
the plain equivalent is documented alongside it, so nothing about nicloj
depends on that library:

```sh
javac -d target/classes $(find nicloj -name '*.java')
```

Only `org.clojure/clojure` is a runtime dependency. Nothing in the source
assumes how `clojure` or `javac` are reached. `dev/in-dev.sh` is a convenience
for the machine this was developed on; it is gitignored and nothing in the
project depends on it. `dev/check.sh` is portable and is tracked.

## Data representations

| Concept | Clojure side | Java side |
| --- | --- | --- |
| Header | `nicloj.header.NiftiHeader` (mutable, but only ever mutated through copying accessors) | same |
| Affine | vector of four vectors of doubles | `double[][]` |
| Voxel array | `nicloj.array.NdArray` | same |
| Voxel values | `double`, always | `double[]` |
| Datatype | keyword, e.g. `:int16` | `DataType` enum |
| Enumerated codes | keyword, e.g. `:aligned`, `:mm` | `int` plus a `Codes` lookup table |
| Orientation | vector of `[out-axis flip]` pairs, `nil` for a dropped axis | `double[][]` with `NaN` for a dropped axis |

`nicloj.core.linalg` and `nicloj.core.ndarray` are thin bridges that convert
between the two columns; they exist so no other namespace has to.

### Mutability

`NiftiHeader` is a mutable struct because that is what a binary layout wants.
Every Clojure accessor that changes a field copies first, so the API is
functional from the outside:

```clojure
(def h (nii/new-header :shape [4 5 6]))
(nii/zooms (nii/set-zooms h [2 2 2]))  ;=> [2.0 2.0 2.0]
(nii/zooms h)                          ;=> [1.0 1.0 1.0]
```

### Voxel arrays are column-major

`NdArray` stores its `double[]` with the first axis varying fastest, matching
the NIfTI file itself, so decoding is a straight copy. `nd/values` yields
elements in that on-disk order; `nd/nested` gives the row-major nesting you
would write by hand. Everything is doubles: an `int16` image costs four times
its file size in memory, which is the price of one uniform representation.

### Laziness

`load` reads the header and returns immediately; the voxel block sits behind a
`delay` that reopens the file on first use. `raw-data` and `fdata` force it,
`loaded?` reports whether it has happened, and `:eager? true` forces it up
front. This matters for the common case of reading a few thousand headers to
find the volumes you want.

### Scaling

An image loaded from disk keeps its voxels exactly as stored, with
`scl_slope`/`scl_inter` in the header; `fdata` applies them, `raw-data` does
not. An image built in memory from data has those fields reset to `1`/`0`,
because the data supplied is already in real-world units. That invariant is why
`nii/image` and `nii/with-data` call `clear-scaling`.
