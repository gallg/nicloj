(ns nicloj.array-test
  "Voxel array layout, reshaping and the byte codec."
  (:require [clojure.test :refer [deftest is testing]]
            [nicloj.core.ndarray :as nd])
  (:import (java.nio ByteOrder)
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

(deftest large-64-bit-values-are-approximate
  (testing "voxels are doubles, so int64/uint64 magnitudes past 2^53 are lossy"
    (let [raw (.array (doto (java.nio.ByteBuffer/allocate 8) (.order le) (.putLong -2)))
          decoded (Codec/decode raw (int-array [1]) DataType/UINT64 le 1.0 0.0)
          reencoded (Codec/encode decoded DataType/UINT64 le 1.0 0.0)]
      (is (= [(Math/pow 2 64)] (nd/values decoded))
          "2^64-2 is not a double, so it decodes as 2^64")
      (is (not= (seq raw) (seq reencoded))
          "and re-encoding cannot recover the original bits"))))

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
