(ns nicloj.header-test
  "Header defaults, field accessors and the binary layout."
  (:require [clojure.test :refer [deftest is testing]]
            [nicloj.core.header :as hdr]
            [nicloj.core.linalg :as la]
            [nicloj.fixtures :refer [all-close?]])
  (:import (nicloj.header Extension NiftiError NiftiHeader)))

(def ^:private affine
  [[2.0 0.0 0.0 -10.0] [0.0 2.0 0.0 -20.0] [0.0 0.0 2.5 -15.0] [0.0 0.0 0.0 1.0]])

(deftest defaults-match-the-nifti-reference
  (let [h (hdr/new-header)]
    (is (= 1 (hdr/version h)))
    (is (hdr/single-file? h))
    (is (not (hdr/big-endian? h)))
    (is (= :float32 (hdr/data-dtype h)))
    (is (= [] (hdr/shape h)) "dim[0] starts at zero")
    (is (= [1.0 0.0] (hdr/slope-inter h)))
    (is (= :unknown (hdr/qform-code h)))
    (is (= :unknown (hdr/sform-code h)))
    (is (= "n+1" (.magic h)))
    (is (= 348 (alength (.toBytes h)))))
  (is (= 540 (alength (.toBytes (hdr/new-header :version 2)))))
  (is (= "ni1" (.magic (hdr/new-header :single? false)))))

(deftest accessors-are-pure
  (let [h (hdr/new-header :shape [4 5 6])
        changed (hdr/set-zooms h [1.5 1.5 3.0])]
    (is (= [1.0 1.0 1.0] (hdr/zooms h)) "the original is untouched")
    (is (= [1.5 1.5 3.0] (hdr/zooms changed)))))

(deftest shape-and-zooms
  (let [h (hdr/new-header :shape [4 5 6 7] :zooms [2 2 2 0.8])]
    (is (= [4 5 6 7] (hdr/shape h)))
    (is (all-close? [2.0 2.0 2.0 0.8] (hdr/zooms h)))
    (is (= [8 8 8] (hdr/shape (hdr/set-shape h [8 8 8])))))
  (testing "zooms must match the number of dimensions"
    (is (thrown? NiftiError (hdr/set-zooms (hdr/new-header :shape [4 5 6]) [1 1])))
    (is (thrown? NiftiError (hdr/set-zooms (hdr/new-header :shape [4]) [-1]))))
  (is (thrown? NiftiError (hdr/new-header :shape (repeat 8 2))) "at most seven dimensions"))

(deftest datatype-names
  (is (= :int16 (hdr/dtype :int16)))
  (is (= :int16 (hdr/dtype "i2")))
  (is (= :int16 (hdr/dtype 4)))
  (is (= :uint8 (hdr/data-dtype (hdr/set-data-dtype (hdr/new-header) :mask))))
  (is (thrown? NiftiError (hdr/dtype :complex256))))

(deftest sform-and-qform-round-trip
  (let [h (hdr/set-affine (hdr/new-header :shape [4 5 6]) affine)]
    (is (= :aligned (hdr/sform-code h)))
    (is (= :unknown (hdr/qform-code h)) "nicloj writes qform values but leaves the code unset")
    (is (all-close? affine (hdr/sform h)))
    (is (all-close? affine (hdr/qform h)))
    (is (all-close? affine (hdr/best-affine h)) "sform wins when both are present"))
  (testing "explicit codes"
    (let [h (hdr/set-qform (hdr/new-header :shape [4 5 6]) affine :scanner)]
      (is (= :scanner (hdr/qform-code h)))
      (is (all-close? affine (hdr/best-affine h)) "qform is used when sform is unset")))
  (testing "clearing a form"
    (let [h (hdr/set-sform (hdr/set-affine (hdr/new-header :shape [4 5 6]) affine) nil)]
      (is (= :unknown (hdr/sform-code h))))))

(deftest base-affine-is-the-last-resort
  (let [h (hdr/new-header :shape [4 5 6] :zooms [2 2 2])]
    (is (= :unknown (hdr/qform-code h)))
    (is (= :unknown (hdr/sform-code h)))
    (is (all-close? (hdr/base-affine h) (hdr/best-affine h)))
    (is (all-close? [-2.0 2.0 2.0] (mapv #(nth (nth (hdr/base-affine h) %) %) [0 1 2]))
        "the Analyze left-right flip is applied")))

(deftest base-affine-handles-a-dimensionless-header
  (let [h (hdr/new-header)]
    (is (= [] (hdr/shape h)))
    (is (= [1.0] (hdr/zooms h)) "zooms reports a single 1.0 when dim[0] is 0")
    (is (all-close? [[-1.0 0.0 0.0 0.0] [0.0 1.0 0.0 0.0]
                     [0.0 0.0 1.0 0.0] [0.0 0.0 0.0 1.0]]
                    (hdr/base-affine h)))
    (is (all-close? (hdr/base-affine h) (hdr/best-affine h)))))

(deftest scaling-fields
  (let [h (hdr/set-slope-inter (hdr/new-header) 0.25 -3.5)]
    (is (= [0.25 -3.5] (hdr/slope-inter h)))
    (is (= [1.0 0.0] (hdr/slope-inter (hdr/clear-scaling h)))))
  (is (thrown? NiftiError (hdr/set-slope-inter (hdr/new-header) 0.0 0.0)))
  (testing "a zero or NaN slope means no scaling"
    (let [h (hdr/copy (hdr/new-header))]
      (set! (.sclSlope h) 0.0)
      (is (nil? (hdr/slope-inter h))))))

(deftest units-and-dim-info
  (let [h (hdr/set-xyzt-units (hdr/new-header) :mm :sec)]
    (is (= [:mm :sec] (hdr/xyzt-units h)))
    (is (= [:micron :msec] (hdr/xyzt-units (hdr/set-xyzt-units h :micron :msec))))
    (is (= [:unknown :unknown] (hdr/xyzt-units (hdr/set-xyzt-units h nil nil)))))
  (let [h (hdr/set-dim-info (hdr/new-header) 0 1 2)]
    (is (= [0 1 2] (hdr/dim-info h)))
    (is (= [nil 2 nil] (hdr/dim-info (hdr/set-dim-info h nil 2 nil))))
    (is (thrown? NiftiError (hdr/set-dim-info h 3 nil nil)))))

;; --------------------------------------------------------------- binary layout

(defn- reread [^NiftiHeader h]
  (hdr/from-bytes (hdr/to-bytes h)))

(deftest binary-round-trip-preserves-every-field
  (doseq [version [1 2]]
    (testing (str "NIfTI-" version)
      (let [h (-> (hdr/new-header :version version :shape [4 5 6 7] :dtype :int16)
                  (hdr/set-zooms [1.5 1.5 3.0 2.0])
                  (hdr/set-affine affine)
                  (hdr/set-slope-inter 0.25 -3.5)
                  (hdr/set-xyzt-units :mm :msec)
                  (hdr/set-dim-info 0 1 2)
                  (hdr/set-descrip "round trip"))
            back (reread h)]
        (is (= (dissoc (hdr/->map h) :vox-offset)
               (dissoc (hdr/->map back) :vox-offset))))))
  (testing "big-endian headers survive"
    (let [h (hdr/copy (hdr/new-header :shape [3 4 5]))]
      (set! (.order h) java.nio.ByteOrder/BIG_ENDIAN)
      (let [back (reread (hdr/set-affine h affine))]
        (is (hdr/big-endian? back))
        (is (all-close? affine (hdr/sform back)))))))

(deftest string-fields-are-truncated-not-corrupted
  (let [long-text (apply str (repeat 200 "x"))
        back (reread (hdr/set-descrip (hdr/new-header) long-text))]
    (is (= 80 (count (hdr/descrip back))))
    (is (= (subs long-text 0 80) (hdr/descrip back)))))

(deftest extensions-survive-a-round-trip
  (let [h (hdr/copy (hdr/new-header :shape [2 2 2]))
        payload (.getBytes "nicloj extension payload")]
    (.add (.extensions h) (Extension. 6 payload))
    (is (= 384 (.minVoxOffset h)) "348 header + 4 extender + 32 padded block")
    (let [back (reread h)]
      (is (= 1 (count (hdr/extensions back))))
      (is (= 6 (.code (first (hdr/extensions back)))))
      (is (= (seq payload) (seq (.content (first (hdr/extensions back)))))))))

(deftest bad-headers-are-rejected
  (is (thrown? NiftiError (hdr/from-bytes (byte-array 400))) "sizeof_hdr is not 348 or 540")
  (is (thrown? NiftiError (hdr/from-bytes (byte-array 2))) "too short")
  (testing "a plausible size but no magic"
    (let [buf (.toBytes (hdr/new-header))]
      (java.util.Arrays/fill buf 344 348 (byte 0))
      (is (thrown? NiftiError (hdr/from-bytes buf))))))

(deftest nifti1-rejects-dimensions-it-cannot-hold
  (let [h (hdr/new-header :shape [40000 2 2])]
    (is (thrown? NiftiError (hdr/to-bytes h)) "dim is int16 in NIfTI-1")
    (is (= [40000 2 2] (hdr/shape (reread (hdr/set-shape (hdr/new-header :version 2)
                                                         [40000 2 2]))))
        "but NIfTI-2 has room")))

(deftest describe-mentions-the-essentials
  (let [text (hdr/describe (hdr/set-affine (hdr/new-header :shape [4 5 6]) affine))]
    (is (re-find #"NIfTI-1" text))
    (is (re-find #"\[4 5 6\]" text))
    (is (re-find #"float32" text))))
