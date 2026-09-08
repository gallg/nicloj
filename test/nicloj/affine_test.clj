(ns nicloj.affine-test
  "Matrix algebra, quaternions and orientation bookkeeping."
  (:require [clojure.test :refer [deftest is testing]]
            [nicloj.core.linalg :as la]
            [nicloj.core.ndarray :as nd]
            [nicloj.core.orientation :as ornt]
            [nicloj.fixtures :refer [all-close? close?]])
  (:import (nicloj.affine Affines Mat Quaternions Svd)
           (nicloj.header NiftiError)))

(def ^:private diag-affine
  [[2.0 0.0 0.0 -10.0] [0.0 3.0 0.0 -20.0] [0.0 0.0 4.0 -30.0] [0.0 0.0 0.0 1.0]])

(defn- rotation
  "Rotation about the z axis by `deg`, as a 3x3 matrix."
  [deg]
  (let [a (Math/toRadians deg) c (Math/cos a) s (Math/sin a)]
    [[c (- s) 0.0] [s c 0.0] [0.0 0.0 1.0]]))

(def ^:private oblique
  (la/from-mat-vec (la/mmul (rotation 25) [[1.5 0 0] [0 1.5 0] [0 0 3.0]])
                   [-31.5 12.25 -8.0]))

;; ----------------------------------------------------------------- matrix ops

(deftest inverse-and-determinant
  (is (all-close? (la/identity 4) (la/mmul diag-affine (la/minverse diag-affine))))
  (is (close? 24.0 (la/mdet diag-affine)))
  (is (close? 1.0 (la/mdet (la/identity 4))))
  (is (all-close? (la/identity 4) (la/mmul oblique (la/minverse oblique))))
  (is (thrown? ArithmeticException (la/minverse [[1.0 2.0] [2.0 4.0]]))))

(deftest affine-application
  (is (all-close? [-8.0 -17.0 -26.0] (la/apply-affine diag-affine [1 1 1])))
  (is (all-close? [[-10.0 -20.0 -30.0] [-8.0 -17.0 -26.0]]
                  (la/apply-affine diag-affine [[0 0 0] [1 1 1]])))
  (testing "the inverse affine undoes it"
    (is (all-close? [3.0 4.0 5.0]
                    (la/apply-affine (la/minverse diag-affine)
                                     (la/apply-affine diag-affine [3 4 5]))))))

(deftest mat-vec-decomposition-round-trips
  (doseq [a [diag-affine oblique]]
    (let [[linear translation] (la/to-mat-vec a)]
      (is (= 3 (count linear)))
      (is (all-close? a (la/from-mat-vec linear translation))))))

(deftest voxel-sizes-are-column-norms
  (is (all-close? [2.0 3.0 4.0] (la/voxel-sizes diag-affine)))
  (is (all-close? [1.5 1.5 3.0] (la/voxel-sizes oblique))))

;; -------------------------------------------------------------------- svd

(deftest svd-reconstructs-its-input
  (doseq [m [(rotation 37) [[3.0 1.0 0.0] [1.0 3.0 0.0] [0.0 0.0 2.0]]
             [[1.0 2.0 3.0] [4.0 5.0 6.0] [7.0 8.0 9.0]]]]
    (let [svd (Svd/of (la/arr m))
          s (vec (.s svd))
          reconstructed (la/mmul (la/mmul (la/rows (.u svd)) (la/rows (Mat/diag (double-array s))))
                                 (la/rows (Mat/transpose (.v svd))))]
      (is (all-close? m reconstructed 1e-10) "u * diag(s) * v^T")
      (is (= s (vec (sort > s))) "singular values are descending"))))

(deftest polar-of-a-rotation-is-the-rotation
  (let [r (rotation 40)]
    (is (all-close? r (la/rows (.polar (Svd/of (la/arr r)) 0.0)) 1e-12)))
  (testing "shear is projected away"
    (let [sheared [[1.0 0.3 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]]
          p (la/rows (.polar (Svd/of (la/arr sheared)) 0.0))]
      (is (all-close? (la/identity 3) (la/mmul p (la/rows (Mat/transpose (la/arr p)))) 1e-12)
          "the result is orthogonal"))))

;; ------------------------------------------------------------- quaternions

(deftest quaternion-matrix-round-trip
  (doseq [deg [0 15 90 179]]
    (let [r (rotation deg)
          q (Quaternions/fromMatrix (la/arr r))]
      (is (>= (aget q 0) 0.0) "the real part is kept non-negative")
      (is (all-close? r (la/rows (Quaternions/toMatrix q)) 1e-12)))))

(deftest fill-positive-recovers-the-real-part
  (is (all-close? [1.0 0.0 0.0 0.0] (vec (Quaternions/fillPositive (double-array [0 0 0]) 1e-7))))
  (is (all-close? [0.0 1.0 0.0 0.0] (vec (Quaternions/fillPositive (double-array [1 0 0]) 1e-7))))
  (is (thrown? NiftiError (Quaternions/fillPositive (double-array [1 1 1]) 1e-7))))

(deftest qform-round-trip
  (doseq [affine [diag-affine oblique]]
    (let [q (Affines/toQform (la/arr affine))
          back (Affines/fromQform (.quaternB q) (.quaternC q) (.quaternD q)
                                  (.qoffset q) (.zooms q) (.qfac q))]
      (is (all-close? affine (la/rows back) 1e-10))
      (is (not (.sheared q))))))

(deftest qform-handles-left-handed-affines
  (let [flipped (assoc-in diag-affine [0 0] -2.0)
        q (Affines/toQform (la/arr flipped))]
    (is (= -1.0 (.qfac q)) "handedness moves into pixdim[0]")
    (is (all-close? [2.0 3.0 4.0] (vec (.zooms q))))
    (is (all-close? flipped
                    (la/rows (Affines/fromQform (.quaternB q) (.quaternC q) (.quaternD q)
                                                (.qoffset q) (.zooms q) (.qfac q)))
                    1e-10))))

(deftest qform-reports-dropped-shear
  (let [sheared (assoc-in diag-affine [0 1] 0.9)]
    (is (.sheared (Affines/toQform (la/arr sheared))))))

(deftest shape-zoom-affine-centres-the-image
  (let [aff (la/rows (Affines/shapeZoomAffine (int-array [4 5 6]) (double-array [2 2 2]) true))]
    (is (all-close? [[-2.0 0.0 0.0 3.0] [0.0 2.0 0.0 -4.0] [0.0 0.0 2.0 -5.0] [0.0 0.0 0.0 1.0]]
                    aff))
    (is (all-close? [0.0 0.0 0.0] (la/apply-affine aff [1.5 2.0 2.5]))
        "the image centre maps to the origin")))

;; ------------------------------------------------------------- orientations

(deftest io-orientation-of-simple-affines
  (is (= [[0 1] [1 1] [2 1]] (ornt/io-orientation (la/identity 4))))
  (is (= ["R" "A" "S"] (ornt/axcodes (la/identity 4))))
  (is (= ["L" "A" "S"] (ornt/axcodes (assoc-in (la/identity 4) [0 0] -1.0))))
  (testing "a transposed affine reports the permutation"
    (let [swapped [[0.0 2.0 0.0 0.0] [1.0 0.0 0.0 0.0] [0.0 0.0 1.0 0.0] [0.0 0.0 0.0 1.0]]]
      (is (= [[1 1] [0 1] [2 1]] (ornt/io-orientation swapped)))
      (is (= ["A" "R" "S"] (ornt/axcodes swapped))))))

(deftest axcodes-round-trip
  (doseq [codes [["R" "A" "S"] ["L" "P" "I"] ["P" "S" "R"] ["A" "L" "S"]]]
    (is (= codes (ornt/ornt->axcodes (ornt/axcodes->ornt codes)))))
  (is (thrown? NiftiError (ornt/axcodes->ornt ["R" "A" "Q"]))))

(deftest ornt-transform-composes
  (let [lia (ornt/axcodes->ornt ["L" "I" "A"])
        ras (ornt/axcodes->ornt ["R" "A" "S"])]
    (is (= ras (ornt/ornt-transform ras ras)))
    (testing "transforming then applying gives the target codes"
      (let [t (ornt/ornt-transform lia ras)]
        (is (= 3 (count t)))
        (is (every? #(#{1 -1} (second %)) t))))))

(deftest inv-ornt-aff-inverts-flips
  (let [shape [4 5 6]
        flip-x [[0 -1] [1 1] [2 1]]
        aff (ornt/inv-ornt-aff flip-x shape)]
    (is (all-close? [[-1.0 0.0 0.0 3.0] [0.0 1.0 0.0 0.0] [0.0 0.0 1.0 0.0] [0.0 0.0 0.0 1.0]]
                    aff))
    (is (all-close? (la/identity 4) (la/mmul aff aff) 1e-12) "flipping twice is a no-op"))
  (is (thrown? NiftiError (ornt/inv-ornt-aff [[nil nil] [1 1] [2 1]] [1 1 1]))))

(deftest apply-orientation-moves-voxels-to-match
  (let [arr (nd/array [2 3 4] (range 24))]
    (testing "the canonical orientation leaves data alone"
      (is (= (nd/values arr) (nd/values (ornt/apply-orientation arr ornt/canonical)))))
    (testing "a flip on axis 0 reverses that axis"
      (is (= (nd/values (nd/flip arr 0))
             (nd/values (ornt/apply-orientation arr [[0 -1] [1 1] [2 1]])))))
    (testing "a permutation transposes"
      (is (= (nd/values (nd/transpose arr [1 0 2]))
             (nd/values (ornt/apply-orientation arr [[1 1] [0 1] [2 1]])))))
    (testing "extra axes are carried along"
      (let [arr4 (nd/reshape arr [2 3 2 2])]
        (is (= [3 2 2 2] (nd/shape (ornt/apply-orientation arr4 [[1 1] [0 1] [2 1]]))))))))
