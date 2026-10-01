(ns nicloj.array-test
  "Voxel array layout, reshaping and the byte codec."
  (:require [clojure.test :refer [deftest is testing]]
            [nicloj.core.ndarray :as nd])
  (:import (java.nio ByteOrder)
           (nicloj.core.ndarray ArrayView)
           (nicloj.array Codec)
           (nicloj.header DataType NiftiError)))

(def ^:private a234 (nd/array [2 3 4] (range 24)))
(def ^:private le ByteOrder/LITTLE_ENDIAN)

(deftest column-major-layout
  (testing "the first axis varies fastest, as on disk"
    (is (= 0.0 (nd/value a234 [0 0 0])))
    (is (= 1.0 (nd/value a234 [1 0 0])))
    (is (= 2.0 (nd/value a234 [0 1 0])))
    (is (= 6.0 (nd/value a234 [0 0 1])))
    (is (= 23.0 (nd/value a234 [1 2 3]))))
  (is (= (map double (range 24)) (nd/values a234)) "values follows the flat order"))

(deftest nested-conversion-is-row-major
  (let [a (nd/nested->array [[1 2 3] [4 5 6]])]
    (is (= [2 3] (nd/shape a)))
    (is (= 2.0 (nd/value a [0 1])))
    (is (= 4.0 (nd/value a [1 0])))
    (is (= [[1.0 2.0 3.0] [4.0 5.0 6.0]] (nd/nested a))))
  (is (= (nd/nested a234) (nd/nested (nd/coerce (nd/nested a234))))
      "nested and coerce are inverses"))

(deftest reshape-and-squeeze
  (is (= [6 4] (nd/shape (nd/reshape a234 [6 4]))))
  (is (= (nd/values a234) (nd/values (nd/reshape a234 [24]))))
  (is (= [2 3 4] (nd/shape (nd/squeeze (nd/reshape a234 [2 1 3 4 1])))))
  (is (thrown? IllegalArgumentException (nd/reshape a234 [5 5]))))

(deftest transpose-permutes-axes
  (let [t (nd/transpose a234 [2 0 1])]
    (is (= [4 2 3] (nd/shape t)))
    (is (= (nd/value a234 [1 2 3]) (nd/value t [3 1 2])))
    (is (= (nd/nested a234) (nd/nested (nd/transpose t [1 2 0]))) "and is invertible")))

(deftest flip-reverses-one-axis
  (let [f (nd/flip a234 1)]
    (is (= [2 3 4] (nd/shape f)))
    (is (= (nd/value a234 [1 0 2]) (nd/value f [1 2 2])))
    (is (= (nd/nested a234) (nd/nested (nd/flip f 1))) "and is its own inverse")))

(deftest slice-extracts-strided-blocks
  (let [s (nd/slice a234 [nil [1 3 1] [0 4 2]])]
    (is (= [2 2 2] (nd/shape s)))
    (is (= (nd/value a234 [0 1 0]) (nd/value s [0 0 0])))
    (is (= (nd/value a234 [1 2 2]) (nd/value s [1 1 1]))))
  (is (thrown? IndexOutOfBoundsException (nd/slice a234 [[0 9 1] nil nil]))))

(deftest concat-joins-along-an-axis
  (let [j (nd/concat [a234 a234] 0)]
    (is (= [4 3 4] (nd/shape j)))
    (is (= (nd/value a234 [1 2 3]) (nd/value j [1 2 3])))
    (is (= (nd/value a234 [1 2 3]) (nd/value j [3 2 3]))))
  (is (= [2 3 8] (nd/shape (nd/concat [a234 a234] 2))))
  (is (thrown? IllegalArgumentException
               (nd/concat [a234 (nd/array [3 3 4] (range 36))] 2))))

(deftest finite-range-ignores-nan-and-infinity
  (is (= [0.0 23.0] (nd/finite-range a234)))
  (is (= [1.0 3.0] (nd/finite-range (nd/array [4] [1 ##NaN 3 ##Inf]))))
  (is (nil? (nd/finite-range (nd/array [2] [##NaN ##NaN])))))

;; --------------------------------------------------------------------- codec

(def ^:private real-types
  [:uint8 :int8 :int16 :uint16 :int32 :uint32 :int64 :uint64 :float32 :float64])

(defn- probe-values
  "Extremes plus a few ordinary values, all exactly representable in `dt`."
  [^DataType dt]
  (if (.isInteger dt)
    (filterv #(<= (.minValue dt) % (.maxValue dt))
             [(.minValue dt) -1.0 0.0 1.0 42.0 (.maxValue dt)])
    [-1024.0 -1.0 0.0 0.25 1.0 1024.0]))

(deftest every-real-type-round-trips-in-both-byte-orders
  (doseq [label real-types]
    (testing (name label)
      (let [dt (DataType/fromLabel (name label))
            values (probe-values dt)
            arr (nd/array [(count values)] values)]
        (doseq [order [ByteOrder/LITTLE_ENDIAN ByteOrder/BIG_ENDIAN]]
          (let [bytes (Codec/encode arr dt order 1.0 0.0)
                back (Codec/decode bytes (int-array [(count values)]) dt order 1.0 0.0)]
            (is (= (* (count values) (.itemSize dt)) (alength bytes)) "byte count")
            (is (nd/close? arr back 0.0) (str order " round trip"))))))))

(deftest integer-encoding-rounds-half-to-even-and-clamps
  (let [arr (nd/array [5] [-40000.0 -0.5 0.5 1.5 40000.0])
        back (Codec/decode (Codec/encode arr DataType/INT16 le 1.0 0.0)
                           (int-array [5]) DataType/INT16 le 1.0 0.0)]
    (is (= [-32768.0 0.0 0.0 2.0 32767.0] (nd/values back)))))

(deftest nan-encodes-as-zero
  (let [arr (nd/array [2] [##NaN 5.0])
        back (Codec/decode (Codec/encode arr DataType/INT16 le 1.0 0.0)
                           (int-array [2]) DataType/INT16 le 1.0 0.0)]
    (is (= [0.0 5.0] (nd/values back)))))

(deftest scaling-is-applied-and-inverted
  (let [arr (nd/array [4] [-3.5 -3.25 0.0 10.0])
        bytes (Codec/encode arr DataType/INT16 le 0.25 -3.5)
        back (Codec/decode bytes (int-array [4]) DataType/INT16 le 0.25 -3.5)]
    (is (nd/close? arr back 1e-12))))

(deftest auto-scale-fits-the-target-range
  (testing "integral data already in range needs no scaling"
    (is (= [1.0 0.0] (vec (Codec/autoScale (nd/array [3] [-5 0 300]) DataType/INT16 true)))))
  (testing "float targets are never scaled"
    (is (= [1.0 0.0] (vec (Codec/autoScale (nd/array [2] [1e30 -1e30])
                                           DataType/FLOAT32 true)))))
  (testing "float data is mapped onto the integer type range"
    (let [arr (nd/array [3] [-1.0 0.25 1.0])
          [slope inter] (vec (Codec/autoScale arr DataType/INT16 true))
          back (Codec/decode (Codec/encode arr DataType/INT16 le slope inter)
                             (int-array [3]) DataType/INT16 le slope inter)]
      (is (not= 1.0 slope) "scaling was needed")
      (is (nd/close? arr back 1e-4))))
  (testing "out-of-range integers are rescaled to fit"
    (let [arr (nd/array [2] [0.0 1e6])
          [slope inter] (vec (Codec/autoScale arr DataType/UINT8 true))
          back (Codec/decode (Codec/encode arr DataType/UINT8 le slope inter)
                             (int-array [2]) DataType/UINT8 le slope inter)]
      (is (nd/close? arr back 1e4)))))

(deftest large-64-bit-values-survive-a-round-trip
  (testing "voxels keep their on-disk type, so int64/uint64 past 2^53 are exact"
    (let [raw (.array (doto (java.nio.ByteBuffer/allocate 8) (.order le) (.putLong -2)))
          decoded (Codec/decode raw (int-array [1]) DataType/UINT64 le 1.0 0.0)
          reencoded (Codec/encode decoded DataType/UINT64 le 1.0 0.0)]
      (is (= :uint64 (nd/dtype decoded)) "the array holds uint64, not doubles")
      (is (= (seq raw) (seq reencoded))
          "re-encoding recovers the original bits exactly")
      (is (= [(Math/pow 2 64)] (nd/values decoded))
          "reading one out as a double still rounds, since 2^64-2 is not a double")))
  (testing "the exact path needs both sides to agree on the type"
    (let [raw (.array (doto (java.nio.ByteBuffer/allocate 8) (.order le) (.putLong -2)))
          decoded (Codec/decode raw (int-array [1]) DataType/UINT64 le 1.0 0.0)]
      (is (not= (seq raw) (seq (Codec/encode decoded DataType/INT64 le 1.0 0.0)))
          "widening through a double to a different type is still lossy"))))

(deftest datatype-lookup
  (is (= DataType/INT16 (DataType/fromCode 4)))
  (is (= DataType/FLOAT64 (DataType/fromLabel "f8")))
  (is (= DataType/UINT8 (DataType/fromLabel "MASK")))
  (is (= 32 (.bitpix DataType/FLOAT32)))
  (is (thrown? NiftiError (DataType/fromCode 1536)) "float128 is not supported")
  (is (thrown? NiftiError (DataType/fromLabel "banana"))))

(deftest unsupported-voxel-types-are-rejected
  (is (thrown? NiftiError (Codec/encode (nd/array [1] [0]) DataType/RGB24 le 1.0 0.0)))
  (is (thrown? NiftiError
               (Codec/decode (byte-array 3) (int-array [1]) DataType/RGB24 le 1.0 0.0))))

(deftest truncated-input-is-rejected
  (is (thrown? NiftiError (Codec/decode (byte-array 4) (int-array [10]) DataType/INT16 le 1.0 0.0))))

(deftest printing-abbreviates-large-arrays
  (let [big (nd/array [50 50] (range 2500))]
    (testing "a big array shows its ends, not every element"
      (let [s (str big)]
        (is (re-find #"\.\.\." s) "the elided middle is marked")
        (is (< (count s) (/ (nd/size big) 4))
            "the repr stays far smaller than one line per element")
        (is (= (count (re-seq #"\n" s))
               (count (re-seq #"\n" (str (nd/array [500 500] (range 250000))))))
            "the repr's height is fixed, not proportional to the array")))
    (testing "views print the same way and stay cheap"
      (is (re-find #"\.\.\." (str (nd/values big))))
      (is (re-find #"\.\.\." (str (nd/nested big))))))
  (testing "a small array prints in full"
    (let [s (str (nd/array [2 3] (range 6)))]
      (is (not (re-find #"\.\.\." s)))
      (is (every? #(re-find (re-pattern (str % "\\.")) s) (range 6))))))

(deftest views-are-lazy-but-behave-as-data
  (let [^ArrayView v (nd/values a234)
        ^ArrayView n (nd/nested a234)]
    (testing "printing does not realise the elements"
      (is (string? (str v)))
      (is (string? (str n)))
      (is (not (realized? (.-d v))))
      (is (not (realized? (.-d n)))))
    (testing "and the data is there when asked for"
      (is (= (nd/size a234) (count v)))
      (is (= (nd/value a234 [1 2 3]) (nth v (dec (count v))) (get-in n [1 2 3])))
      (is (== (reduce + (range 24)) (reduce + v)))
      (is (= (vec (nd/values a234)) (into [] v)))
      (is (realized? (.-d v))))))

(deftest voxels-keep-their-on-disk-type
  (testing "decoding yields the file's type, not doubles"
    (doseq [[dt kw] [[DataType/UINT8 :uint8] [DataType/INT8 :int8]
                     [DataType/INT16 :int16] [DataType/UINT16 :uint16]
                     [DataType/INT32 :int32] [DataType/UINT32 :uint32]
                     [DataType/INT64 :int64] [DataType/UINT64 :uint64]
                     [DataType/FLOAT32 :float32] [DataType/FLOAT64 :float64]]]
      (let [^DataType dt dt
            raw (byte-array (* 4 (.itemSize dt)))
            arr (Codec/decode raw (int-array [4]) dt le 1.0 0.0)]
        (is (= kw (nd/dtype arr)) (str "decode keeps " kw))
        (is (= [0.0 0.0 0.0 0.0] (nd/values arr)) (str "and reads as doubles for " kw)))))
  (testing "structural operations carry the type through"
    (let [arr (Codec/decode (byte-array 24) (int-array [2 3 2]) DataType/INT16 le 1.0 0.0)]
      (doseq [[label a] [["reshape" (nd/reshape arr [12])]
                         ["squeeze" (nd/squeeze (nd/reshape arr [2 3 2 1]))]
                         ["transpose" (nd/transpose arr [2 1 0])]
                         ["flip" (nd/flip arr 0)]
                         ["slice" (nd/slice arr [[0 1 1] nil nil])]
                         ["concat" (nd/concat [arr arr] 0)]]]
        (is (= :int16 (nd/dtype a)) (str label " keeps int16")))))
  (testing "scaling widens, because scaled values are no longer integers"
    (let [arr (Codec/decode (byte-array 8) (int-array [4]) DataType/INT16 le 1.0 0.0)]
      (is (= :int16 (nd/dtype (nd/scaled arr 1.0 0.0))) "identity scaling is a no-op")
      (is (= :float64 (nd/dtype (nd/scaled arr 2.0 1.0)))))))

(deftest flat-access-mutates-in-place
  (testing "flat/setFlat walk the buffer in column-major order without boxing an index"
    (let [orig (nd/array [2 3] (range 6))
          a (nd/array [2 3] (range 6))]
      (is (= (nd/values a) (mapv #(.flat a %) (range (nd/size a)))))
      (dotimes [i (nd/size a)] (.setFlat a i (- (.flat a i))))
      (is (= (mapv - (nd/values orig)) (nd/values a)))
      (is (= (- (nd/value orig [1 2])) (nd/value a [1 2]))
          "and address the same elements as the index accessors")))
  (testing "the buffer keeps its type, so writes are narrowed to it"
    (let [a (Codec/decode (byte-array 4) (int-array [2]) DataType/INT16 le 1.0 0.0)]
      (.setFlat a 0 3.7)
      (is (= :int16 (nd/dtype a)))
      (is (= 3.0 (nd/value a [0])) "3.7 stored into int16 truncates"))))

(deftest bad-permutations-and-ragged-input-are-rejected
  (is (thrown? IllegalArgumentException (nd/transpose a234 [0 0 1])))
  (is (thrown? IllegalArgumentException (nd/transpose a234 [0 1 3])))
  (is (thrown? IllegalArgumentException (nd/coerce [[1 2] [3]])))
  (is (thrown? IllegalArgumentException (nd/coerce [[1 2] [3 4 5]])))
  (is (thrown? IllegalArgumentException (nd/coerce [[1 2] 3]))))

(deftest views-grow-and-hash-like-the-vectors-they-equal
  (let [v (nd/values (nd/array [3] [1 2 3]))]
    (is (= [1.0 2.0 3.0 4.0] (conj v 4.0)))
    (is (= [1.0 2.0 3.0 4.0 5.0] (into v [4.0 5.0])))
    (is (= (hash [1.0 2.0 3.0]) (hash v)))
    (is (contains? #{[1.0 2.0 3.0]} v)))
  (is (= 3.5 (nd/nested (nd/array [] [3.5]))) "0-d nests as its element"))

(deftest close-treats-equal-infinities-as-equal
  (let [x (nd/array [2] [##Inf ##-Inf])]
    (is (nd/close? x x))
    (is (not (nd/close? x (nd/array [2] [##-Inf ##Inf]))))))

(deftest bad-arguments-get-clear-errors
  (is (thrown-with-msg? IllegalArgumentException #"axis 3 out of range" (nd/flip a234 3)))
  (is (thrown-with-msg? IllegalArgumentException #"axis 3 out of range" (nd/concat [a234 a234] 3)))
  (is (thrown-with-msg? IllegalArgumentException #"must be integers" (nd/value a234 [1.5 0 0])))
  (is (thrown? IllegalArgumentException (nd/coerce [[1 nil]])))
  (is (thrown? IllegalArgumentException (nd/coerce [[1 "x"]])))
  (is (= [1 3 4] (nd/shape (nd/slice a234 [[1 nil] nil nil]))) "nil inside a triple"))
