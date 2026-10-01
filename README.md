<h1><img src="icon.png" alt="" height="64" align="middle"> nicloj</h1>

NIfTI-1 and NIfTI-2 reading and writing in Clojure and Java, written from
scratch and tested against [nibabel](https://nipy.org/nibabel/).

Handles `.nii`, `.nii.gz` and NIfTI `.hdr`/`.img` pairs (not the older Analyze
7.5 format), either byte order, every integer datatype plus `float32` and
`float64`, qform/sform affines, orientations and header extensions.

## Use

```clojure
(require '[nicloj.api.image :as nii])

(def img (nii/load "brain.nii.gz"))

(nii/shape img)                 ;=> [64 64 30 120]
(nii/zooms img)                 ;=> [3.0 3.0 4.0 2.5]
(nii/affine img)                ;=> [[-3.0 0.0 0.0 96.0] ...]
(nii/axcodes (nii/affine img))  ;=> ["L" "A" "S"]
(nii/voxel img [32 32 15 0])    ;=> 412.5
(nii/fdata img)                 ; the whole volume, read on first use

(-> img
    nii/as-closest-canonical
    (nii/slice-image [[10 54] [10 54] nil nil])
    (nii/save "cropped.nii.gz" :dtype :int16))
```

`nicloj.api.image` re-exports the whole API; the modules underneath
(`nicloj.core.nifti`, `nicloj.core.header`, ...) can be required directly too.

## Build and test

Needs `clojure` and `javac` on `PATH`. Compile the Java modules once:

```sh
clojure -T:build javac
```

To test the implementation against nibabel, run `./dev/check.sh`. It needs a
Python environment with nibabel and nilearn.

## Layout

```
nicloj/header/   Java   header structs, datatypes, extensions
nicloj/array/    Java   voxel array and byte codec
nicloj/affine/   Java   matrices, quaternions, orientations
nicloj/io/       Clj    file names and gzip streams
nicloj/core/     Clj    header API, image type, reader/writer, ops
nicloj/api/      Clj    the public facade
```

`progress/` has the design notes, current state and remaining work.

## License

MIT, see [LICENSE](LICENSE). The NIfTI layout comes from the public-domain
`nifti1.h`/`nifti2.h` reference headers; some algorithms are reimplemented from
nibabel (MIT), attributed in the source.
