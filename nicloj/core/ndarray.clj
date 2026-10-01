(ns nicloj.core.ndarray
  "Clojure view of `nicloj.array.NdArray`.

  Voxel arrays are column-major: the first axis varies fastest, as in the NIfTI
  file itself. `values` therefore yields elements in on-disk order, while
  `nested` gives the row-major nesting you would write by hand.

  Elements always read as doubles, but are stored at their on-disk width --
  see `dtype`. Array operations carry that type through."
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
              (let [depth (count idx)]
                (cond
                  (and (sequential? x) (< depth (count shape)) (= (count x) (nth shape depth)))
                  (doseq [[i child] (map-indexed vector x)] (walk child (conj idx i)))
                  (and (number? x) (= depth (count shape)))
                  (.set a (double x) (int-array idx))
                  :else (throw (IllegalArgumentException.
                                (str "nested data must be numbers in equal-length rows; bad entry at "
                                     idx " for shape " shape))))))]
      (walk nested [])
      a)))

(defn coerce
  "Coerce `data` to an `NdArray`: passes one through, otherwise reads nested
  sequences row-major."
  ^NdArray [data]
  (if (instance? NdArray data) data (nested->array data)))

(defn dtype
  "The voxel type the elements are held in, as a keyword."
  [^NdArray a]
  (keyword (.label (.dtype a))))

(defn shape [^NdArray a] (vec (.shape a)))
(defn ndim [^NdArray a] (.ndim a))
(defn size [^NdArray a] (.size a))

(defn value
  "The voxel at `idx`, a sequence of per-axis indices."
  ^double [^NdArray a idx]
  (when-not (every? integer? idx)
    (throw (IllegalArgumentException. (str "indices must be integers, got " idx))))
  (.get a (int-array idx)))

(deftype ArrayView [^NdArray array ^String header ^clojure.lang.Delay d]
  clojure.lang.Sequential
  clojure.lang.IPersistentCollection
  (seq [_] (seq @d))
  (count [_] (count @d))
  (cons [_ x] (conj @d x))
  (empty [_] [])
  (equiv [_ o] (= @d o))
  clojure.lang.IHashEq
  (hasheq [_] (hash @d))
  clojure.lang.Indexed
  (nth [_ i] (nth @d i))
  (nth [_ i not-found] (nth @d i not-found))
  clojure.lang.ILookup
  (valAt [_ k] (get @d k))
  (valAt [_ k not-found] (get @d k not-found))
  clojure.lang.IReduceInit
  (reduce [_ f init] (reduce f init @d))
  Iterable
  (iterator [_] (clojure.lang.RT/iter @d))
  Object
  (hashCode [_] (.hashCode ^Object @d))
  (equals [_ o] (.equals ^Object @d o))
  (toString [_] (str header "\n" (.layout array))))

(defmethod print-method ArrayView [v ^java.io.Writer w]
  (.write w (str v)))

(defn values
  "All elements in column-major (on-disk) order.

  The elements are realised on first use, and printing shows only the ends of
  the array, so handing a large one back to a REPL is cheap."
  [^NdArray a]
  (->ArrayView (.reshape a (int-array [(.size a)]))
               (str "#nicloj/values[" (.size a) "]")
               (delay (vec (.toDoubleArray a)))))

(defn nested
  "Elements as row-major nested vectors, realised on first use. A 0-d array
  gives its one element, as numpy's `tolist` does."
  [^NdArray a]
  (if (zero? (.ndim a))
    (.get a (int-array 0))
    (->ArrayView a (str "#nicloj/nested" (shape a))
                 (delay
                   (let [shape (shape a)]
                     (letfn [(walk [idx depth]
                               (if (= depth (count shape))
                                 (.get a (int-array idx))
                                 (mapv #(walk (conj idx %) (inc depth)) (range (nth shape depth)))))]
                       (walk [] 0)))))))

(defn reshape ^NdArray [^NdArray a shape] (.reshape a (int-array shape)))
(defn squeeze ^NdArray [^NdArray a] (.squeeze a))
(defn transpose ^NdArray [^NdArray a perm] (.transpose a (int-array perm)))
(defn flip ^NdArray [^NdArray a axis] (.flip a (int axis)))

(defn slice
  "Extract the sub-block given by one `[start stop step]` triple per axis;
  `nil` in place of a triple keeps the whole axis."
  ^NdArray [^NdArray a specs]
  (let [dims (shape a)
        specs (map-indexed (fn [i [start stop step]]
                             [(or start 0) (or stop (nth dims i)) (or step 1)])
                           specs)]
    (.slice a
            (int-array (map first specs))
            (int-array (map second specs))
            (int-array (map last specs)))))

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

(defmethod print-method NdArray [^NdArray a ^java.io.Writer w]
  (.write w (.toString a)))
