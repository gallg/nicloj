(ns nicloj.core.linalg
  "Clojure view of the `nicloj.affine` matrix code.

  Matrices are nested vectors of doubles here and `double[][]` inside Java;
  `rows` and `arr` convert between the two."
  (:refer-clojure :exclude [identity])
  (:import (nicloj.affine Affines Mat)))

(defn arr
  "Convert a nested sequence of numbers into a `double[][]`."
  ^"[[D" [m]
  (if (instance? (Class/forName "[[D") m)
    m
    (into-array (Class/forName "[D") (map double-array m))))

(defn rows
  "Convert a `double[][]` into a vector of vectors of doubles."
  [^"[[D" m]
  (mapv vec m))

(defn identity
  "The `n` by `n` identity matrix."
  [n]
  (rows (Mat/identity n)))

(defn mmul
  "Matrix product of `a` and `b`."
  [a b]
  (rows (Mat/mul (arr a) (arr b))))

(defn minverse
  "Matrix inverse of the square matrix `a`."
  [a]
  (rows (Mat/inverse (arr a))))

(defn mdet
  "Determinant of the square matrix `a`."
  ^double [a]
  (Mat/det (arr a)))

(defn close?
  "True when `a` and `b` agree elementwise to within `atol` (default 1e-8)."
  ([a b] (close? a b 1e-8))
  ([a b atol] (Mat/close (arr a) (arr b) atol)))

(defn apply-affine
  "Map `point` (or each row of a sequence of points) through the affine `a`."
  [a point]
  (if (number? (first point))
    (vec (Affines/apply (arr a) (double-array point)))
    (mapv vec (Affines/applyAll (arr a) (arr point)))))

(defn voxel-sizes
  "Voxel edge lengths implied by an affine: the norms of its rotation columns."
  [a]
  (vec (Affines/voxelSizes (arr a))))

(defn from-mat-vec
  "Assemble an affine from a linear block and a translation vector."
  [linear translation]
  (rows (Affines/fromMatVec (arr linear) (double-array translation))))

(defn to-mat-vec
  "Split an affine into `[linear translation]`, the inverse of `from-mat-vec`."
  [a]
  (let [m (mapv vec a)
        n (dec (count m))]
    [(mapv #(subvec % 0 (dec (count %))) (subvec m 0 n))
     (mapv peek (subvec m 0 n))]))
