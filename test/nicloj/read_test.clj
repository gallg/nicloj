(ns nicloj.read-test
  "Check nicloj's reader against nibabel, field by field, over the whole corpus."
  (:require [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [nicloj.api.image :as nii]
            [nicloj.core.header :as hdr]
            [nicloj.fixtures :as fix :refer [all-close? close?]])
  (:import (nicloj.array NdArray)))

(def ^:private stat-keys [:min :max :sum :mean])

(defn- stats-of
  "Min, max, sum and mean over finite voxels, the way the manifest records them.
  Sums compensated so they stay comparable with numpy's pairwise summation."
  [^NdArray arr]
  (let [^doubles d (.data arr)
        n (alength d)]
    (loop [i 0, mn Double/POSITIVE_INFINITY, mx Double/NEGATIVE_INFINITY
           sum 0.0, err 0.0, count 0]
      (if (= i n)
        {:min mn :max mx :sum sum :mean (/ sum count)}
        (let [v (aget d i)]
          (if (Double/isFinite v)
            (let [y (- v err)
                  t (+ sum y)]
              (recur (inc i) (min mn v) (max mx v) t (- (- t sum) y) (inc count)))
            (recur (inc i) mn mx sum err count)))))))

(defn- check-stats [want arr label]
  (is (all-close? (mapv want stat-keys) (mapv (stats-of arr) stat-keys) 1e-9) label))

(defn- check-samples [samples arr label]
  (doseq [[idx value] samples]
    (is (close? value (nii/value arr idx) 1e-12) (str label " at " idx))))

(deftest header-fields-match-nibabel
  (doseq [{:keys [file] :as want} @fix/manifest]
    (testing file
      (let [img (nii/load (fix/corpus-file file))
            h (nii/header img)
            m (nii/header->map h)]
        (is (= (:version want) (nii/version h)) "NIfTI version")
        (is (= (:byte-order want) (:byte-order m)) "byte order")
        (is (= (:single-file? want) (nii/single-file? h)) "single file vs pair")
        (is (= (:shape want) (nii/shape img)) "shape")
        (is (= (:dtype want) (nii/data-dtype img)) "datatype")
        (is (all-close? (:zooms want) (nii/zooms img)) "zooms")
        (is (= (:vox-offset want) (nii/data-offset h)) "vox_offset")
        (is (all-close? (:scl want) [(:scl-slope m) (:scl-inter m)]) "scl_slope/scl_inter")
        (is (= (:qform-code want) (nii/qform-code h)) "qform_code")
        (is (= (:sform-code want) (nii/sform-code h)) "sform_code")
        (is (= (:xyzt-units want) (nii/xyzt-units h)) "xyzt_units")
        (is (= (:dim-info want) (nii/dim-info h)) "dim_info")
        (is (= (:descrip want) (nii/descrip h)) "descrip")
        (is (= (:intent-code want) (first (nii/intent h))) "intent_code")
        (is (= (:slice-code want) (nii/slice-code h)) "slice_code")
        (is (= (:n-extensions want) (count (nii/extensions h))) "extension count")))))

(deftest affines-match-nibabel
  (doseq [{:keys [file] :as want} @fix/manifest]
    (testing file
      (let [img (nii/load (fix/corpus-file file))
            h (nii/header img)]
        (is (all-close? (:affine want) (nii/affine img)) "best affine")
        (is (all-close? (:qform want) (nii/qform h)) "qform")
        (is (all-close? (:sform want) (nii/sform h)) "sform")
        (is (= (:axcodes want) (nii/axcodes (nii/affine img))) "axis codes")
        (is (= (:ornt want) (nii/io-orientation (nii/affine img))) "io orientation")))))

(deftest voxels-match-nibabel
  (doseq [{:keys [file] :as want} @fix/manifest]
    (testing file
      (let [data (nii/fdata (nii/load (fix/corpus-file file)))]
        (check-stats (:stats want) data "voxel statistics")
        (check-samples (:samples want) data "voxel")))))

(deftest canonical-reorientation-matches-nibabel
  (doseq [{:keys [file canonical]} @fix/manifest]
    (testing file
      (let [canon (nii/as-closest-canonical (nii/load (fix/corpus-file file)))
            data (nii/fdata canon)]
        (is (= (:shape canonical) (nii/shape canon)) "canonical shape")
        (is (all-close? (:affine canonical) (nii/affine canon)) "canonical affine")
        (is (= (:axcodes canonical) (nii/axcodes (nii/affine canon))) "canonical axcodes")
        (check-stats (:stats canonical) data "canonical statistics")
        (check-samples (:samples canonical) data "canonical voxel")))))

(deftest voxels-are-read-lazily
  (let [img (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))]
    (is (not (nii/loaded? img)) "no voxels read yet")
    (is (= [7 8 9] (nii/shape img)) "shape comes from the header alone")
    (is (some? (nii/fdata img)))
    (is (nii/loaded? img) "voxels cached once read"))
  (is (nii/loaded? (nii/load (fix/corpus-file "nifti1-3d-float32.nii") :eager? true))))

(deftest reads-header-alone
  (let [h (nii/read-header (fix/corpus-file "nifti2-3d-float32.nii"))]
    (is (= 2 (nii/version h)))
    (is (= [6 6 6] (hdr/shape h)))
    (is (= :float32 (hdr/data-dtype h)))))

(deftest either-half-of-a-pair-resolves
  (let [from-hdr (nii/load (fix/corpus-file "nifti1-pair.hdr"))
        from-img (nii/load (fix/corpus-file "nifti1-pair.img"))]
    (is (= (nii/shape from-hdr) (nii/shape from-img)))
    (is (nii/array-close? (nii/fdata from-hdr) (nii/fdata from-img))))
  (let [gz (nii/load (fix/corpus-file "nifti1-pair-gz.img.gz"))]
    (is (= [4 5 6] (nii/shape gz)))))

(deftest pair-halves-may-disagree-about-gzip
  (let [img (nii/load (fix/corpus-file "nifti1-pair.hdr"))]
    (nii/save img (fix/out-file "mixed.hdr") :scaling :keep)
    (nii/save img (fix/out-file "mixed-gz.hdr.gz") :scaling :keep)
    ;; Leave a plain .hdr beside a gzipped .img and nothing else.
    (io/copy (io/file (fix/out-file "mixed-gz.img.gz")) (io/file (fix/out-file "mixed.img.gz")))
    (io/delete-file (fix/out-file "mixed.img"))
    (testing "the gzipped image is found from the plain header name"
      (is (nii/array-close? (nii/fdata img)
                            (nii/fdata (nii/load (fix/out-file "mixed.hdr"))) 0.0)))
    (testing "the plain header is found from the gzipped image name"
      (is (nii/array-close? (nii/fdata img)
                            (nii/fdata (nii/load (fix/out-file "mixed.img.gz"))) 0.0)))))

(deftest rejects-unreadable-input
  (testing "unknown extension"
    (is (thrown? clojure.lang.ExceptionInfo (nii/load "somewhere/image.mgz"))))
  (testing "not a NIfTI header"
    (let [path (fix/out-file "garbage.nii")]
      (spit path (apply str (repeat 400 "x")))
      (is (thrown? nicloj.header.NiftiError (nii/load path))))))
