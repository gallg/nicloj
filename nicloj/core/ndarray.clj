(ns nicloj.core.ndarray
  "Clojure view of `nicloj.array.NdArray`.

  Voxel arrays are column-major: the first axis varies fastest, as in the NIfTI
  file itself. `values` therefore yields elements in on-disk order, while
  `nested` gives the row-major nesting you would write by hand."
  (:refer-clojure :exclude [concat])
  (:import (nicloj.array NdArray)))

(defn array
  "Build an array of `shape` from `flat`, a column-major sequence of numbers."
  ^NdArray [shape flat]
  (NdArray. (int-array shape) (double-array flat)))

(defn zeros ^NdArray [shape] (NdArray. (int-array shape)))

(defn- nested-shape [x]
  (if (sequential? x) (into [(count x)] (nested-shape (first x))) []))

(defn nested->array
  "Build an array from row-major nested sequences, e.g. `[[1 2] [3 4]]`."
  ^NdArray [nested]
  (let [shape (nested-shape nested)
        a (NdArray. (int-array shape))]
    (letfn [(walk [x idx]
              (if (sequential? x)
                (doseq [[i child] (map-indexed vector x)] (walk child (conj idx i)))
                (.set a (double x) (int-array idx))))]
      (walk nested [])
      a)))

(defn coerce
  "Coerce `data` to an `NdArray`: passes one through, otherwise reads nested
  sequences row-major."
  ^NdArray [data]
  (if (instance? NdArray data) data (nested->array data)))

(defn shape [^NdArray a] (vec (.shape a)))
(defn ndim [^NdArray a] (.ndim a))
(defn size [^NdArray a] (.size a))

(defn value
  "The voxel at `idx`, a sequence of per-axis indices."
  ^double [^NdArray a idx]
  (.get a (int-array idx)))

(defn values
  "All elements in column-major (on-disk) order."
  [^NdArray a]
  (vec (.data a)))

(defn nested
  "Elements as row-major nested vectors."
  [^NdArray a]
  (let [shape (shape a)]
    (letfn [(walk [idx depth]
              (if (= depth (count shape))
                (.get a (int-array idx))
                (mapv #(walk (conj idx %) (inc depth)) (range (nth shape depth)))))]
      (walk [] 0))))

(defn reshape ^NdArray [^NdArray a shape] (.reshape a (int-array shape)))
(defn squeeze ^NdArray [^NdArray a] (.squeeze a))
(defn transpose ^NdArray [^NdArray a perm] (.transpose a (int-array perm)))
(defn flip ^NdArray [^NdArray a axis] (.flip a (int axis)))

(defn slice
  "Extract the sub-block given by one `[start stop step]` triple per axis;
  `nil` in place of a triple keeps the whole axis."
  ^NdArray [^NdArray a specs]
  (let [dims (shape a)
        specs (map-indexed (fn [i s] (or s [0 (nth dims i) 1])) specs)]
    (.slice a
            (int-array (map first specs))
            (int-array (map second specs))
            (int-array (map #(nth % 2 1) specs)))))

(defn concat
  "Join arrays end to end along `axis`."
  ^NdArray [arrays axis]
  (NdArray/concat (into-array NdArray arrays) (int axis)))

(defn scaled
  "Elementwise `v * slope + inter`."
  ^NdArray [^NdArray a slope inter]
  (.scaled a (double slope) (double inter)))

(defn finite-range
  "`[min max]` over finite elements, or nil when there are none."
  [^NdArray a]
  (let [[mn mx] (vec (.finiteRange a))]
    (when (<= mn mx) [mn mx])))

(defn close?
  "True when `a` and `b` have the same shape and agree to within `atol`."
  ([a b] (close? a b 1e-9))
  ([^NdArray a ^NdArray b atol] (.closeTo a b (double atol))))
