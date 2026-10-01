(ns nicloj.ops-test
  "Image construction, coordinate mapping and the whole-image operations."
  (:require [clojure.test :refer [deftest is testing]]
            [nicloj.api.image :as nii]
            [nicloj.fixtures :as fix :refer [all-close?]]))

(def ^:private affine
  [[2.0 0.0 0.0 -10.0] [0.0 2.0 0.0 -20.0] [0.0 0.0 2.5 -15.0] [0.0 0.0 0.0 1.0]])

(defn- volume [shape]
  (nii/image (nii/array shape (range (apply * shape))) affine))

;; --------------------------------------------------------------- construction

(deftest building-an-image
  (let [img (volume [2 3 4])]
    (is (nii/image? img))
    (is (= [2 3 4] (nii/shape img)))
    (is (= 3 (nii/ndim img)))
    (is (= :float32 (nii/data-dtype img)) "the default on-disk type")
    (is (all-close? [2.0 2.0 2.5] (nii/zooms img)) "zooms come from the affine")
    (is (all-close? affine (nii/affine img)))
    (is (= 23.0 (nii/voxel img [1 2 3])))
    (is (= [1.0 0.0] (nii/slope-inter (nii/header img)))
        "data given directly is already in real-world units")))

(deftest building-from-nested-vectors
  (let [img (nii/image [[[1 2] [3 4]] [[5 6] [7 8]]] (nii/eye 4))]
    (is (= [2 2 2] (nii/shape img)))
    (is (= 1.0 (nii/voxel img [0 0 0])))
    (is (= 8.0 (nii/voxel img [1 1 1])))))

(deftest replacing-parts-of-an-image
  (let [img (volume [2 3 4])]
    (testing "with-data updates dim"
      (let [smaller (nii/with-data img (nii/array [2 2] [1 2 3 4]))]
        (is (= [2 2] (nii/shape smaller)))
        (is (all-close? affine (nii/affine smaller)) "the affine is left alone")))
    (testing "with-affine rewrites both forms"
      (let [moved (nii/with-affine img (nii/eye 4))]
        (is (all-close? (nii/eye 4) (nii/affine moved)))
        (is (all-close? (nii/eye 4) (nii/sform (nii/header moved))))
        (is (all-close? (nii/eye 4) (nii/qform (nii/header moved))))))
    (testing "with-dtype only changes what gets written"
      (is (= :int16 (nii/data-dtype (nii/with-dtype img :int16))))
      (is (= 23.0 (nii/voxel (nii/with-dtype img :int16) [1 2 3]))))))

(deftest with-header-preserves-voxel-values
  (let [img (nii/load (fix/corpus-file "nifti1-3d-int16-scaled.nii.gz"))
        value (nii/voxel img [0 0 0])
        scl (nii/slope-inter (nii/header img))]
    (is (= (:scl (fix/entry "nifti1-3d-int16-scaled.nii.gz")) scl))
    (testing "an unrelated header does not silently rescale the stored data"
      (let [swapped (nii/with-header img (nii/new-header :shape [6 7 5] :dtype :int16))]
        (is (= value (nii/voxel swapped [0 0 0])))
        (is (= scl (nii/slope-inter (nii/header swapped))))))
    (testing "other fields of the new header are still adopted"
      (let [tagged (nii/with-header img (nii/set-descrip (nii/new-header :shape [6 7 5]) "tag"))]
        (is (= "tag" (nii/descrip (nii/header tagged))))
        (is (= value (nii/voxel tagged [0 0 0])))))))

;; ---------------------------------------------------------------- coordinates

(deftest voxel-and-world-coordinates
  (let [img (volume [4 5 6])]
    (is (all-close? [-10.0 -20.0 -15.0] (nii/voxel->world img [0 0 0])))
    (is (all-close? [-8.0 -18.0 -12.5] (nii/voxel->world img [1 1 1])))
    (is (all-close? [1.0 2.0 3.0] (nii/world->voxel img (nii/voxel->world img [1 2 3]))))))

;; ----------------------------------------------------------------- operations

(deftest squeeze-drops-trailing-singleton-axes
  (is (= [2 3 4] (nii/shape (nii/squeeze-image (volume [2 3 4 1 1])))))
  (is (= [2 3 4 2] (nii/shape (nii/squeeze-image (volume [2 3 4 2 1])))))
  (testing "the first three axes are never dropped"
    (is (= [1 1 1] (nii/shape (nii/squeeze-image (volume [1 1 1]))))))
  (testing "an image with nothing to squeeze comes back unchanged"
    (let [img (volume [2 3 4])]
      (is (identical? img (nii/squeeze-image img))))))

(deftest concat-images-stacks-and-joins
  (let [a (volume [2 3 4])
        b (volume [2 3 4])]
    (testing "with no axis, a new trailing axis appears"
      (let [joined (nii/concat-images [a b])]
        (is (= [2 3 4 2] (nii/shape joined)))
        (is (= (nii/voxel a [1 2 3]) (nii/voxel joined [1 2 3 0])))
        (is (= (nii/voxel b [1 2 3]) (nii/voxel joined [1 2 3 1])))))
    (testing "with an axis, the images are joined along it"
      (is (= [4 3 4] (nii/shape (nii/concat-images [a b] :axis 0)))))
    (testing "affines must agree"
      (is (thrown? nicloj.header.NiftiError
                   (nii/concat-images [a (nii/with-affine b (nii/eye 4))])))
      (is (some? (nii/concat-images [a (nii/with-affine b (nii/eye 4))]
                                    :check-affines? false))))
    (is (thrown? nicloj.header.NiftiError (nii/concat-images [])))))

(deftest four-to-three-splits-volumes
  (let [img (volume [2 3 4 5])
        parts (nii/four-to-three img)]
    (is (= 5 (count parts)))
    (is (every? #(= [2 3 4] (nii/shape %)) parts))
    (is (= (nii/voxel img [1 2 3 4]) (nii/voxel (last parts) [1 2 3])))
    (is (all-close? affine (nii/affine (first parts))))
    (testing "and rejoining reproduces the original"
      (is (nii/array-close? (nii/fdata img) (nii/fdata (nii/concat-images parts)) 0.0))))
  (is (thrown? nicloj.header.NiftiError (nii/four-to-three (volume [2 3 4])))))

(deftest slice-image-crops-and-moves-the-affine
  (let [img (volume [6 6 6])
        cropped (nii/slice-image img [[2 5] [1 4] nil])]
    (is (= [3 3 6] (nii/shape cropped)))
    (is (= (nii/voxel img [2 1 0]) (nii/voxel cropped [0 0 0])))
    (testing "the same voxel keeps its world position"
      (is (all-close? (nii/voxel->world img [2 1 0]) (nii/voxel->world cropped [0 0 0])))
      (is (all-close? (nii/voxel->world img [4 3 5]) (nii/voxel->world cropped [2 2 5])))))
  (testing "subsampling scales the affine"
    (let [img (volume [6 6 6])
          coarse (nii/slice-image img [[0 6 2] [0 6 2] [0 6 2]])]
      (is (= [3 3 3] (nii/shape coarse)))
      (is (all-close? [4.0 4.0 5.0] (nii/voxel-sizes (nii/affine coarse))))
      (is (all-close? (nii/voxel->world img [2 2 2]) (nii/voxel->world coarse [1 1 1]))))))

;; --------------------------------------------------------------- reorientation

(deftest reorienting-preserves-world-positions
  (doseq [file ["ornt-lia.nii.gz" "ornt-psr.nii.gz" "ornt-als.nii.gz" "oblique.nii.gz"]]
    (testing file
      (let [img (nii/load (fix/corpus-file file))
            canon (nii/as-closest-canonical img)]
        (is (= ["R" "A" "S"] (nii/axcodes (nii/affine canon))))
        (testing "every sampled voxel keeps its world coordinate and value"
          (doseq [idx [[0 0 0] [1 2 3] [2 1 4]]]
            (let [world (nii/voxel->world img idx)
                  back (mapv #(Math/round (double %)) (nii/world->voxel canon world))]
              (is (all-close? world (nii/voxel->world canon back) 1e-6))
              (is (fix/close? (nii/voxel img idx) (nii/voxel canon back) 1e-6)))))))))

(deftest reorienting-is-a-no-op-when-already-canonical
  (let [img (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))]
    (is (identical? img (nii/as-closest-canonical img)))
    (is (identical? img (nii/as-reoriented img nii/canonical-ornt)))))

(deftest reorienting-updates-dim-info
  (let [img (nii/load (fix/corpus-file "nifti1-3d-float32.nii"))
        rotated (nii/as-reoriented img [[1 1] [0 1] [2 1]])]
    (is (= (:dim-info (fix/entry "nifti1-3d-float32.nii")) (nii/dim-info (nii/header img))))
    (is (= (mapv {0 1, 1 0, 2 2} (nii/dim-info (nii/header img)))
           (nii/dim-info (nii/header rotated)))
        "the frequency and phase axes swap with the data")))

(deftest concat-images-names-a-mismatched-shape
  (let [a (nii/image [[[1 2]]] affine)
        b (nii/image [[[1 2 3]]] affine)]
    (is (thrown-with-msg? nicloj.header.NiftiError #"image 1 has shape" (nii/concat-images [a b])))))

(deftest slicing-matches-nibabels-slicer
  (doseq [{:keys [file sliced]} @fix/manifest :when sliced]
    (testing file
      (let [want (nii/load (fix/corpus-file file))
            got (nii/slice-image (nii/load (fix/corpus-file (:source sliced))) (:specs sliced))]
        (is (= (nii/shape want) (nii/shape got)))
        (is (all-close? (nii/affine want) (nii/affine got) 1e-6))
        (is (all-close? (nii/zooms want) (nii/zooms got)))
        (is (nii/array-close? (nii/fdata want) (nii/fdata got) 0.0))))))

(deftest integer-specs-drop-axes-past-the-third
  (let [img (nii/image (nii/array [2 2 2 3 4] (range 96)) affine)
        img (nii/with-header img (nii/set-zooms (nii/header img) [2 2 2.5 3 7]))]
    (testing "the dropped axis takes its voxel size with it"
      (is (= [2 2 2 4] (nii/shape (nii/slice-image img [nil nil nil 1]))))
      (is (all-close? [2 2 2.5 7] (nii/zooms (nii/slice-image img [nil nil nil 1])))))
    (testing "negative indices count from the end"
      (is (nii/array-close? (nii/fdata (nii/slice-image img [nil nil nil 2 3]))
                            (nii/fdata (nii/slice-image img [nil nil nil -1 -1])) 0.0)))
    (doseq [[specs msg] [[[nil nil 1] #"spatial axis 2"]
                         [[nil nil nil 3] #"out of range"]
                         [[nil nil nil -4] #"out of range"]]]
      (is (thrown-with-msg? nicloj.header.NiftiError msg (nii/slice-image img specs))))))
