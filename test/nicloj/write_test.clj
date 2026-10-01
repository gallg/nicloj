(ns nicloj.write-test
  "Round-trip every corpus image through nicloj's writer, in each layout.

  The files land in test-data/out/ together with written.edn, which
  scripts/verify_roundtrip.py then re-checks with nibabel."
  (:require [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [nicloj.api.image :as nii]
            [nicloj.core.header :as hdr]
            [nicloj.fixtures :as fix :refer [all-close?]]
            [nicloj.io.paths :as paths]))

(def variants
  "Layouts every image is rewritten into. `:keep` scaling rewrites the stored
  voxels untouched, so the reloaded values must match exactly."
  [{:name "as-nii" :suffix ".nii" :opts {:scaling :keep}}
   {:name "as-gz" :suffix ".nii.gz" :opts {:scaling :keep}}
   {:name "as-pair" :suffix ".hdr" :opts {:scaling :keep}}
   {:name "as-nifti2" :suffix "-v2.nii" :opts {:scaling :keep :version 2}}])

(defn- stem [file]
  (reduce (fn [s ext] (if (str/ends-with? s ext) (subs s 0 (- (count s) (count ext))) s))
          file
          [".gz" ".nii" ".hdr" ".img"]))

(defn- compare-images
  "Every property a faithful copy has to preserve."
  [original copy label]
  (is (= (nii/shape original) (nii/shape copy)) (str label " shape"))
  (is (= (nii/data-dtype original) (nii/data-dtype copy)) (str label " datatype"))
  (is (all-close? (nii/zooms original) (nii/zooms copy) 1e-6) (str label " zooms"))
  (is (all-close? (nii/affine original) (nii/affine copy) 1e-6) (str label " affine"))
  (is (= (nii/qform-code (nii/header original)) (nii/qform-code (nii/header copy)))
      (str label " qform_code"))
  (is (= (nii/sform-code (nii/header original)) (nii/sform-code (nii/header copy)))
      (str label " sform_code"))
  (is (= (nii/xyzt-units (nii/header original)) (nii/xyzt-units (nii/header copy)))
      (str label " xyzt_units"))
  (is (= (nii/dim-info (nii/header original)) (nii/dim-info (nii/header copy)))
      (str label " dim_info"))
  (is (= (nii/descrip (nii/header original)) (nii/descrip (nii/header copy)))
      (str label " descrip"))
  (is (= (count (nii/extensions (nii/header original))) (count (nii/extensions (nii/header copy))))
      (str label " extensions"))
  (is (nii/array-close? (nii/fdata original) (nii/fdata copy) 0.0)
      (str label " voxels are identical")))

(deftest round-trips-every-corpus-image
  (let [written (atom [])]
    (doseq [{:keys [file]} @fix/manifest
            {:keys [name suffix opts]} variants]
      (testing (str file " -> " name)
        (let [original (nii/load (fix/corpus-file file))
              out (fix/out-file (str (stem file) "--" name suffix))]
          (apply nii/save original out (mapcat identity opts))
          (swap! written conj {:source file :written (.getName (java.io.File. out))
                               :variant name})
          (compare-images original (nii/load out) name))))
    (spit (fix/out-file "written.edn") (pr-str @written))))

(deftest keep-scaling-preserves-raw-voxels-and-scalers
  (let [original (nii/load (fix/corpus-file "nifti1-3d-int16-scaled.nii.gz"))
        out (fix/out-file "keep-scaling.nii")
        copy (nii/load (nii/save original out :scaling :keep))]
    (is (= (:scl (fix/entry "nifti1-3d-int16-scaled.nii.gz")) (nii/slope-inter (nii/header copy))))
    (is (nii/array-close? (nii/raw-data original) (nii/raw-data copy) 0.0))
    (is (thrown? nicloj.header.NiftiError
                 (nii/save original (fix/out-file "keep-clash.nii")
                           :scaling :keep :dtype :float32))
        ":keep cannot also change the datatype")))

(deftest auto-scaling-fits-a-narrower-datatype
  (let [original (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))
        out (fix/out-file "auto-int16.nii")
        copy (nii/load (nii/save original out :dtype :int16))
        [slope _] (nii/slope-inter (nii/header copy))
        [mn mx] (nii/finite-range (nii/fdata original))
        step (/ (- mx mn) 65535.0)]
    (is (= :int16 (nii/data-dtype copy)))
    (is (not= 1.0 slope) "a slope was derived")
    (is (nii/array-close? (nii/fdata original) (nii/fdata copy) step)
        "values survive to within one quantisation step")))

(deftest explicit-scaling-is-honoured
  (let [original (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))
        out (fix/out-file "explicit-scaling.nii")
        copy (nii/load (nii/save original out :dtype :int16 :scaling [0.01 0.0]))]
    (is (all-close? [0.01 0.0] (nii/slope-inter (nii/header copy)))
        "NIfTI-1 keeps the scalers as float32")
    (is (nii/array-close? (nii/fdata original) (nii/fdata copy) 0.005))))

(deftest vox-offset-and-magic-follow-the-layout
  (let [img (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))]
    (is (= 352 (nii/data-offset (nii/read-header (nii/save img (fix/out-file "off1.nii"))))))
    (is (= 0 (nii/data-offset (nii/read-header (nii/save img (fix/out-file "off1.hdr"))))))
    (is (= 544 (nii/data-offset (nii/read-header (nii/save img (fix/out-file "off2.nii")
                                                           :version 2)))))
    (is (= (:vox-offset (fix/entry "with-extension.nii")) (nii/data-offset
                                                           (nii/read-header (nii/save (nii/load (fix/corpus-file "with-extension.nii"))
                                                                                      (fix/out-file "off-ext.nii") :scaling :keep))))
        "extensions push the voxel offset out")))

(deftest gzip-is-chosen-by-the-file-name
  (let [img (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))
        plain (java.io.File. ^String (nii/save img (fix/out-file "gzip-plain.nii")))
        zipped (java.io.File. ^String (nii/save img (fix/out-file "gzip-zipped.nii.gz")))]
    (is (< (.length zipped) (.length plain)))
    (is (nii/array-close? (nii/fdata img) (nii/fdata (nii/load (.getPath zipped))) 0.0))))

(deftest bytes-round-trip-in-memory
  (let [original (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))
        copy (nii/from-bytes (nii/->bytes original :scaling :keep))]
    (compare-images original copy "->bytes")))

(deftest truncated-bytes-are-rejected
  (let [original (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))
        whole ^bytes (nii/->bytes original :scaling :keep)
        short-buf (java.util.Arrays/copyOf whole (- (alength whole) 100))]
    (is (= (nii/voxel original (mapv dec (nii/shape original)))
           (nii/voxel (nii/from-bytes whole) (mapv dec (nii/shape original))))
        "the full buffer round-trips the last voxel")
    (is (thrown? nicloj.header.NiftiError (nii/from-bytes short-buf))
        "a short buffer must fail, not silently read zeros")))

(deftest big-endian-output-is-readable
  (let [original (nii/load (fix/corpus-file "nifti1-3d-int16-bigendian.nii"))
        out (fix/out-file "bigendian-out.nii")
        copy (nii/load (nii/save original out :scaling :keep))]
    (is (nii/big-endian? (nii/header copy)) "byte order is inherited from the header")))

(deftest writing-a-freshly-built-image
  (let [data (nii/array [3 4 5] (range 60))
        affine [[2.0 0.0 0.0 -3.0] [0.0 2.0 0.0 -4.0] [0.0 0.0 2.0 -5.0] [0.0 0.0 0.0 1.0]]
        h (-> (nii/new-header :shape [3 4 5] :dtype :int16)
              (hdr/set-xyzt-units :mm :sec)
              (hdr/set-descrip "built by nicloj"))
        img (nii/image data affine :header h)
        copy (nii/load (nii/save img (fix/out-file "fresh.nii.gz")))]
    (is (= [3 4 5] (nii/shape copy)))
    (is (= :int16 (nii/data-dtype copy)))
    (is (= [:mm :sec] (nii/xyzt-units (nii/header copy))))
    (is (= "built by nicloj" (nii/descrip (nii/header copy))))
    (is (all-close? affine (nii/affine copy)))
    (is (nii/array-close? data (nii/fdata copy) 0.0))))

(deftest nan-into-a-scaled-integer-reads-back-near-zero
  (let [f (fix/out-file "nan-int16.nii")
        _ (nii/save (nii/image [[[##NaN 1000.5 2000.5]]] (nii/eye 4)) f :dtype :int16)
        img (nii/load f)
        [slope _] (nii/slope-inter (nii/header img))]
    (is (<= (Math/abs (nii/voxel img [0 0 0])) slope))))

(deftest keep-refuses-a-dtype-changed-by-with-dtype
  (let [img (nii/load (fix/corpus-file "nifti1-3d-int16-scaled.nii.gz"))]
    (is (thrown? Exception (nii/save (nii/with-dtype img :uint8)
                                     (fix/out-file "keep-uint8.nii") :scaling :keep)))))

(deftest pair-companions-follow-the-case-of-the-name
  (is (= "UP.IMG" (:image (paths/resolve-image "UP.HDR"))))
  (is (= "x.HDR.GZ" (:header (paths/resolve-image "x.IMG.GZ"))))
  (is (= "x.img.GZ" (:image (paths/resolve-image "x.hdr.GZ")))))

(deftest bad-scaling-and-version-are-refused
  (let [img (nii/image [[[0.5 1.5 2.5]]] (nii/eye 4))
        f (fix/out-file "refused.nii")]
    (doseq [s [[0 0] [##NaN 0] [1 ##NaN] :none [1]]]
      (is (thrown? Exception (nii/save img f :dtype :int16 :scaling s)) (str s)))
    (is (thrown? Exception (nii/save (nii/image [[[1e39 1e39]]] (nii/eye 4)) f :dtype :int16)))
    (doseq [v [0 3]]
      (is (thrown? Exception (nii/save img f :version v)) (str v)))))

(deftest explicit-scalers-are-narrowed-before-encoding-in-nifti-1
  (let [f (fix/out-file "explicit-v1.nii")
        _ (nii/save (nii/image [[[2e8 1e8]]] (nii/eye 4)) f :dtype :int32 :scaling [0.1 0])
        img (nii/load f)
        [slope _] (nii/slope-inter (nii/header img))]
    (is (<= (Math/abs (- 2e8 (nii/voxel img [0 0 0]))) (/ slope 2)))))
