(ns nicloj.processing-test
  "Resampling, checked against nibabel.processing over the corpus and against
  what resampling has to mean on its own."
  (:require [clojure.test :refer [deftest is testing]]
            [nicloj.api.image :as nii]
            [nicloj.fixtures :as fix :refer [all-close?]])
  (:import (nicloj.header NiftiError)))

(defn- resample-like
  "Redo the resampling a manifest entry records nibabel doing."
  [{:keys [source how order voxel-sizes cval target]} want]
  (let [src (nii/load (fix/corpus-file source))]
    (case how
      :to-output (nii/resample-to-output src :order order :voxel-sizes (or voxel-sizes 1.0))
      :from-to (nii/resample-from-to src [(nii/shape want) target]
                                     :order order :cval (or cval 0.0)))))

(deftest resampling-matches-nibabel
  (doseq [{:keys [file resampled]} @fix/manifest :when resampled]
    (testing file
      (let [want (nii/load (fix/corpus-file file))
            got (resample-like resampled want)]
        (is (= (nii/shape want) (nii/shape got)))
        (is (all-close? (nii/affine want) (nii/affine got) 1e-6) "affine, to float32 storage")
        (is (= (nii/data-dtype want) (nii/array-dtype (nii/raw-data got)))
            "voxel type carried through")
        ;; Integer results must agree exactly; float ones to rounding noise.
        (is (nii/array-close? (nii/fdata want) (nii/fdata got) 1e-12))))))

(def ^:private affine
  [[2.0 0.0 0.0 -10.0] [0.0 3.0 0.0 -20.0] [0.0 0.0 4.0 -30.0] [0.0 0.0 0.0 1.0]])

(defn- volume [shape]
  (nii/image (nii/array shape (map #(* 1.5 %) (range (reduce * shape)))) affine))

(deftest resampling-onto-its-own-grid-changes-nothing
  (let [img (volume [4 5 6])]
    (doseq [order [0 1]]
      (let [out (nii/resample-from-to img img :order order)]
        ;; Composing the affines leaves ~1e-15 of residue, but no slice may be lost.
        (is (nii/array-close? (nii/fdata img) (nii/fdata out) 1e-12))
        (is (= (nii/affine img) (nii/affine out)))))))

(deftest a-cropped-grid-gives-the-cropped-voxels
  (let [img (volume [6 6 6])
        crop (nii/slice-image img [[1 5] [2 6 2] nil])]
    (doseq [order [0 1]]
      (is (nii/array-close? (nii/fdata crop)
                            (nii/fdata (nii/resample-from-to img crop :order order))
                            1e-12)))))

(deftest points-outside-take-cval-and-halfway-points-average
  (let [img (volume [3 3 3])
        shifted (assoc-in affine [0 3] -9.0)   ; half a voxel along x
        out (nii/resample-from-to img [[4 3 3] shifted] :cval -1.0)]
    (is (fix/close? (/ (+ (nii/voxel img [0 1 1]) (nii/voxel img [1 1 1])) 2)
                    (nii/voxel out [0 1 1])
                    1e-12))
    (is (= -1.0 (nii/voxel out [3 1 1])) "x = 3.5 voxels is past the last one")))

(deftest a-4d-series-resamples-volume-by-volume-onto-a-3d-grid
  (let [img (volume [4 4 3 2])
        target [[3 3 3] (assoc-in affine [1 3] -18.5)]
        out (nii/resample-from-to img target)]
    (is (= [3 3 3 2] (nii/shape out)))
    (doseq [[t vol] (map-indexed vector (nii/four-to-three img))]
      (is (nii/array-close? (nii/fdata (nii/resample-from-to vol target))
                            (nii/reshape (nii/slice (nii/fdata out) [nil nil nil [t (inc t) 1]])
                                         [3 3 3])
                            0.0)
          (str "volume " t)))))

(deftest to-output-covers-every-source-corner
  (let [img (nii/load (fix/corpus-file "oblique.nii.gz"))
        out (nii/resample-to-output img :voxel-sizes 2.0)
        [sx sy sz] (nii/shape img)
        shape (nii/shape out)]
    (is (= [2.0 2.0 2.0] (nii/zooms out)))
    (is (= ["R" "A" "S"] (nii/axcodes (nii/affine out))))
    (doseq [c (for [i [0 (dec sx)] j [0 (dec sy)] k [0 (dec sz)]] [i j k])
            :let [v (nii/world->voxel out (nii/voxel->world img c))]]
      (is (every? true? (map #(<= -1e-6 %1 (+ (dec %2) 1e-6)) v shape)) (str c)))))

(deftest bad-resampling-requests-are-refused
  (let [img (volume [3 3 3])]
    (is (thrown? IllegalArgumentException (nii/resample-from-to img img :order 3)))
    (is (thrown? NiftiError (nii/resample-from-to (nii/image [[1 2] [3 4]] affine) img)))
    (is (thrown? NiftiError (nii/resample-from-to img [[3 3] affine])))
    (is (thrown? NiftiError (nii/resample-to-output img :voxel-sizes [1 0 1])))
    (is (thrown? NiftiError (nii/resample-to-output (volume [3 3 3 2]))))))
