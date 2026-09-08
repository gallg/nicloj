(ns nicloj.api.image
  "The whole of nicloj behind one require.

      (require '[nicloj.api.image :as nii])

      (def img (nii/load \"brain.nii.gz\"))
      (nii/shape img)                  ;=> [64 64 30 120]
      (nii/affine img)                 ;=> [[-3.0 0.0 0.0 96.0] ...]
      (nii/value (nii/fdata img) [32 32 15 0])

      (-> (nii/image data [[2 0 0 -20] [0 2 0 -30] [0 0 2 -10] [0 0 0 1]])
          (nii/save \"out.nii.gz\" :dtype :int16))

  Each name below is re-exported from the module that implements it, so
  `nicloj.core.nifti`, `nicloj.core.image`, `nicloj.core.header`,
  `nicloj.core.ndarray`, `nicloj.core.orientation`, `nicloj.core.funcs` and
  `nicloj.core.linalg` can also be used directly."
  (:refer-clojure :exclude [load])
  (:require [nicloj.core.funcs :as funcs]
            [nicloj.core.header :as hdr]
            [nicloj.core.image :as image]
            [nicloj.core.linalg :as la]
            [nicloj.core.ndarray :as nd]
            [nicloj.core.nifti :as nifti]
            [nicloj.core.orientation :as ornt]))

(defmacro ^:private reexport
  "Define a var here for each name given, carrying over its docstring and
  argument lists. A `[local qualified]` pair renames on the way through."
  [& specs]
  `(do ~@(for [spec specs
               :let [[local target] (if (vector? spec) spec [(symbol (name spec)) spec])]]
           `(do (def ~local ~target)
                (alter-meta! (var ~local) merge
                             (select-keys (meta (var ~target)) [:doc :arglists]))))))

;; Files
(reexport nifti/load nifti/save nifti/read-header nifti/from-bytes
          [->bytes nifti/->bytes])

;; Images
(reexport image/image image/image? image/header image/affine image/shape image/ndim
          image/zooms image/data-dtype image/fdata image/raw-data image/voxel
          image/loaded? image/with-data image/with-affine image/with-header
          image/with-dtype image/describe)

;; Headers
(reexport hdr/new-header hdr/dtype hdr/version hdr/single-file? hdr/big-endian?
          hdr/data-offset hdr/qform hdr/sform hdr/qform-code hdr/sform-code
          hdr/best-affine hdr/base-affine hdr/set-qform hdr/set-sform hdr/set-affine
          hdr/set-shape hdr/set-zooms hdr/set-data-dtype hdr/slope-inter
          hdr/set-slope-inter hdr/clear-scaling hdr/xyzt-units hdr/set-xyzt-units
          hdr/dim-info hdr/set-dim-info hdr/intent hdr/slice-code hdr/descrip
          hdr/set-descrip hdr/extensions
          [header->map hdr/->map]
          [describe-header hdr/describe])

;; Voxel arrays
(reexport nd/array nd/zeros nd/value nd/values nd/nested nd/reshape nd/squeeze
          nd/transpose nd/flip nd/slice nd/finite-range
          [array-shape nd/shape]
          [array-close? nd/close?])

;; Orientation
(reexport [canonical-ornt ornt/canonical]
          ornt/io-orientation ornt/axcodes ornt/axcodes->ornt ornt/ornt->axcodes
          ornt/ornt-transform ornt/inv-ornt-aff ornt/apply-orientation
          ornt/as-reoriented ornt/as-closest-canonical)

;; Whole-image operations
(reexport funcs/squeeze-image funcs/concat-images funcs/four-to-three
          funcs/slice-image funcs/voxel->world funcs/world->voxel)

;; Affine arithmetic
(reexport la/mmul la/minverse la/mdet la/apply-affine la/voxel-sizes
          la/from-mat-vec la/to-mat-vec
          [eye la/identity])
