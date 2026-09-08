# CLAUDE.md

Context for working on nicloj: a from-scratch NIfTI-1/NIfTI-2 reader and writer
in Clojure and Java, validated against nibabel.

## Reaching the toolchain

`clojure` and `javac` live in a distrobox container on this machine, not on the
host. Run anything Clojure or Java through:

```sh
./dev/in-dev.sh 'clojure -M:test'
```

`dev/in-dev.sh` is a convenience for **this machine only** and is gitignored,
so a fresh clone will not have it. Recreate it with:

```sh
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
exec distrobox enter dev -- bash -lc "cd '$PWD' && $*"
```

nicloj itself must keep working wherever `clojure` and `javac` are on `PATH`.
Never bake distrobox, or any absolute toolchain path, into `deps.edn`,
`build.clj`, `dev/check.sh` or the source.

Python (nibabel 5.4.2, nilearn 0.14.1) runs on the **host**, in the active uv
venv at `.venv/`. Both see the same filesystem, so the Clojure suite and the
Python scripts exchange files through `test-data/`.

## The one thing that will bite you

`clojure -M:test` does **not** compile Java. After editing any `.java` file:

```sh
./dev/in-dev.sh 'clojure -T:build javac'
```

or just run `./dev/check.sh`, which does javac, corpus, tests and the nibabel
verification in order. Two confusing test failures during initial development
were nothing but stale classes in `target/classes`.

## Layout

```
nicloj/header/   Java   NiftiHeader, DataType, Extension, Codes, NiftiError
nicloj/array/    Java   NdArray, Codec
nicloj/affine/   Java   Mat, Svd, Quaternions, Affines, Orientations
nicloj/io/       Clj    nicloj.io.paths, nicloj.io.stream
nicloj/core/     Clj    header, ndarray, linalg, image, nifti, orientation, funcs
nicloj/api/      Clj    nicloj.api.image -- the public facade
```

A module is one folder, one feature, one language. Dependencies run downward
only: `api` → `core` → `io`/`affine`/`array` → `header`. The Java modules know
nothing about Clojure.

The project root is the classpath root (`:paths ["." "target/classes"]`), so
`nicloj/core/nifti.clj` is `nicloj.core.nifti` and
`nicloj/header/NiftiHeader.java` is `nicloj.header.NiftiHeader`.

## Invariants worth not breaking

- **Header accessors are pure.** `NiftiHeader` is mutable, but every `set-*` in
  `nicloj.core.header` copies first via the private `edit` helper. Do not add a
  mutating public accessor.
- **`NdArray` is column-major** — first axis fastest, matching the file. `values`
  returns that order; `nested` returns row-major nesting.
- **Voxels are always `double`.** Consistent, and the reason for the memory
  limitation in `progress/05-status.md`.
- **An image built from data has identity scaling.** `fdata` always applies
  `scl_slope`/`scl_inter`, so data handed to `image` or `with-data` must have
  those fields reset — both call `hdr/clear-scaling`. Skipping that scales twice.
- **`nicloj/affine` is pure maths** and must not depend on `NdArray`.
  `apply-orientation` lives in `nicloj.core.orientation` for that reason.
- **No reflection.** `clojure -M -e "(set! *warn-on-reflection* true) (require
  'nicloj.api.image)"` is clean; keep it that way.
- **Only `org.clojure/clojure` at runtime.** `tools.build` is build-time only,
  and the plain `javac` equivalent is documented so nothing depends on it.

## Testing

Two directions, both must pass:

```sh
./dev/check.sh
```

1. `scripts/gen_testdata.py` builds 29 fixtures with nibabel plus
   `test-data/manifest.edn` recording what nibabel reports for each.
2. `clojure -M:test` (80 tests, ~3800 assertions) asserts nicloj matches the
   manifest field by field, and rewrites every fixture into
   `test-data/out/`.
3. `scripts/verify_roundtrip.py` has nibabel re-read all 116 outputs and compare
   them to their originals — voxels at `atol=0`, plus 25 header fields.

When adding a feature, add its fixture to the generator and its expectation to
the manifest, so both directions cover it. Never hard-code a value copied out
of the corpus into a test — assert the relationship instead, or derive the
number from the data. Two tests did that and broke the moment the corpus
changed.

The corpus is reproducible: `gen_testdata.py` seeds from sha256, so generating
it twice, or in a fresh clone, gives byte-identical files and manifest. Keep it
that way; the earlier `hash()`-based seed made the suite differ run to run and
hid a genuine defect.

**nibabel normalises two things away on load** — it blanks `vox_offset` and
moves `scl_slope`/`scl_inter` into its array proxy. Never use
`hdr.get_slope_inter()` or `hdr['vox_offset']` from a loaded image as an
expectation; the generator reads those from the raw bytes for exactly this
reason.

## Where to read next

`progress/` has the full picture: `01-architecture.md` for the module
boundaries and data representations, `02-nifti-format.md` for the on-disk
layout and every format decision, `03-nibabel-parity.md` for the API mapping
and the deliberate differences, `04-testing.md` for the harness,
`05-status.md` for what works and what does not, `06-todo.md` for what is
next, `07-log.md` for how it got here.

The largest open item is item 1 in `06-todo.md`: keeping voxels in their
on-disk type instead of widening everything to `double`. It changes the
`fdata`/`raw-data` contract, so it is worth doing before other work builds on
the current representation.
