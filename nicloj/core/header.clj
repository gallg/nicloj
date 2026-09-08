(ns nicloj.core.header
  "Reading and updating NIfTI header fields from Clojure.

  Headers are `nicloj.header.NiftiHeader` instances. Accessors are pure; every
  `set-*` returns a modified copy and leaves its argument alone. Enumerated
  fields are keywords (`:int16`, `:aligned`, `:mm`), affines are vectors of
  four vectors of four doubles."
  (:require [clojure.string :as str]
            [nicloj.core.linalg :as la])
  (:import (nicloj.affine Affines)
           (nicloj.header Codes DataType Extension NiftiError NiftiHeader)))

;; ---------------------------------------------------------------- code names

(defn- ->kw [label] (keyword (str/replace label " " "-")))
(defn- ->label [kw] (str/replace (name kw) "-" " "))

(defn- code->kw [table code] (->kw (Codes/label table code)))
(defn- kw->code [table kw] (Codes/code table (->label kw)))

(defn dtype
  "Normalise a datatype given as a keyword, string or NIfTI code to a keyword."
  [t]
  (cond
    (keyword? t) (keyword (.label (DataType/fromLabel (name t))))
    (string? t) (keyword (.label (DataType/fromLabel t)))
    (integer? t) (keyword (.label (DataType/fromCode (int t))))
    (instance? DataType t) (keyword (.label ^DataType t))
    :else (throw (NiftiError. (str "cannot read a datatype from " (pr-str t))))))

(defn- ^DataType data-type [t]
  (if (instance? DataType t) t (DataType/fromLabel (name (dtype t)))))

;; ---------------------------------------------------------------- construction

(defn new-header
  "Build a header. Options:

    :version  1 (default) or 2
    :shape    voxel grid shape
    :dtype    datatype keyword, `:float32` by default
    :zooms    voxel sizes, one per dimension
    :single?  true (default) for a `.nii` image, false for a `.hdr`/`.img` pair"
  ^NiftiHeader [& {:keys [version shape dtype zooms single?]
                   :or {version 1 single? true}}]
  (let [h (NiftiHeader/of version)]
    (set! (.singleFile h) (boolean single?))
    (when shape (.setShape h (int-array shape)))
    (when dtype (.setDataType h (data-type dtype)))
    (when zooms (.setZooms h (double-array zooms)))
    h))

(defn copy ^NiftiHeader [^NiftiHeader h] (.copy h))

(defn- edit
  "Apply the mutating `f` to a copy of `h` and return the copy."
  ^NiftiHeader [^NiftiHeader h f]
  (let [c (.copy h)]
    (f c)
    c))

;; ---------------------------------------------------------------- shape, type

(defn version [^NiftiHeader h] (.version h))
(defn shape [^NiftiHeader h] (vec (.shape h)))
(defn zooms [^NiftiHeader h] (vec (.zooms h)))
(defn data-dtype [^NiftiHeader h] (dtype (.datatype h)))
(defn single-file? [^NiftiHeader h] (.singleFile h))
(defn big-endian? [^NiftiHeader h] (= java.nio.ByteOrder/BIG_ENDIAN (.order h)))
(defn data-offset [^NiftiHeader h] (.voxOffset h))
(defn extensions [^NiftiHeader h] (vec (.extensions h)))

(defn set-shape [^NiftiHeader h shape]
  (edit h #(.setShape ^NiftiHeader % (int-array shape))))

(defn set-zooms [^NiftiHeader h zooms]
  (edit h #(.setZooms ^NiftiHeader % (double-array zooms))))

(defn set-data-dtype [^NiftiHeader h t]
  (edit h #(.setDataType ^NiftiHeader % (data-type t))))

;; -------------------------------------------------------------------- affines

(defn qform-code [^NiftiHeader h] (code->kw Codes/XFORM (.qformCode h)))
(defn sform-code [^NiftiHeader h] (code->kw Codes/XFORM (.sformCode h)))

(defn qform
  "The affine encoded by the quaternion fields, whatever `qform_code` says."
  [^NiftiHeader h]
  (la/rows (Affines/fromQform (.quaternB h) (.quaternC h) (.quaternD h)
                              (double-array [(.qoffsetX h) (.qoffsetY h) (.qoffsetZ h)])
                              (double-array (take 3 (drop 1 (.pixdim h))))
                              (.qfac h))))

(defn sform
  "The affine held in the `srow_*` fields, whatever `sform_code` says."
  [^NiftiHeader h]
  (la/rows (Affines/fromSrow (.srow h))))

(defn base-affine
  "Fallback affine for a header with neither qform nor sform: diagonal, centred
  on the image, with the Analyze left-right flip."
  [^NiftiHeader h]
  (let [shape (.shape h)]
    ;; zooms() reports (1.0) for a dimensionless header, so trim it to match.
    (la/rows (Affines/shapeZoomAffine
              shape (java.util.Arrays/copyOf (.zooms h) (alength shape)) true))))

(defn best-affine
  "The affine nicloj uses for an image: sform if set, else qform, else the
  centred fallback."
  [^NiftiHeader h]
  (cond
    (not= 0 (.sformCode h)) (sform h)
    (not= 0 (.qformCode h)) (qform h)
    :else (base-affine h)))

(defn set-qform
  "Store `affine` in the quaternion fields, also setting `pixdim[0..3]`.

  Any shear is dropped, since the qform can only express rotation, zoom and
  handedness. `code` defaults to `:aligned` when previously unset."
  ([h affine] (set-qform h affine nil))
  ([^NiftiHeader h affine code]
   (edit h (fn [^NiftiHeader c]
             (set! (.qformCode c)
                   (int (cond code (kw->code Codes/XFORM code)
                              (nil? affine) 0
                              (zero? (.qformCode c)) 2
                              :else (.qformCode c))))
             (when affine
               (let [q (Affines/toQform (la/arr affine))
                     pixdim (.pixdim c)]
                 (set! (.quaternB c) (.quaternB q))
                 (set! (.quaternC c) (.quaternC q))
                 (set! (.quaternD c) (.quaternD q))
                 (set! (.qoffsetX c) (aget (.qoffset q) 0))
                 (set! (.qoffsetY c) (aget (.qoffset q) 1))
                 (set! (.qoffsetZ c) (aget (.qoffset q) 2))
                 (aset pixdim 0 (.qfac q))
                 (dotimes [i 3] (aset pixdim (inc i) (aget (.zooms q) i)))))))))

(defn set-sform
  "Store `affine` in the `srow_*` fields. `code` defaults to `:aligned` when
  previously unset."
  ([h affine] (set-sform h affine nil))
  ([^NiftiHeader h affine code]
   (edit h (fn [^NiftiHeader c]
             (set! (.sformCode c)
                   (int (cond code (kw->code Codes/XFORM code)
                              (nil? affine) 0
                              (zero? (.sformCode c)) 2
                              :else (.sformCode c))))
             (when affine
               (let [rows (mapv vec affine)
                     srow (.srow c)]
                 (dotimes [r 3]
                   (let [^doubles row (aget srow r)]
                     (dotimes [k 4]
                       (aset row k (double (get-in rows [r k]))))))))))))

(defn set-affine
  "Put `affine` in both forms the way nicloj writes new images: sform
  `:aligned`, qform values with code `:unknown`."
  [h affine]
  (-> h (set-sform affine :aligned) (set-qform affine :unknown)))

;; ------------------------------------------------------------------- scaling

(defn slope-inter
  "The `[slope inter]` pair to apply to raw voxels, or nil when the header
  asks for no scaling."
  [^NiftiHeader h]
  (let [slope (.sclSlope h)
        inter (.sclInter h)]
    (cond
      (or (zero? slope) (not (Double/isFinite slope))) nil
      (not (Double/isFinite inter))
      (throw (NiftiError. (str "valid scl_slope but invalid scl_inter " inter)))
      :else [slope inter])))

(defn set-slope-inter [^NiftiHeader h slope inter]
  (when (or (zero? slope) (not (Double/isFinite slope)))
    (throw (NiftiError. "scl_slope must be finite and non-zero")))
  (edit h (fn [^NiftiHeader c]
            (set! (.sclSlope c) (double slope))
            (set! (.sclInter c) (double inter)))))

(defn clear-scaling
  "Reset `scl_slope`/`scl_inter` to the identity, as befits voxel data already
  held in real-world units."
  [h]
  (set-slope-inter h 1.0 0.0))

;; ---------------------------------------------------------------- other codes

(defn xyzt-units
  "Spatial and temporal units as `[xyz t]` keywords."
  [^NiftiHeader h]
  (let [[xyz t] (vec (Codes/splitUnits (.xyztUnits h)))]
    [(code->kw Codes/UNITS xyz) (code->kw Codes/UNITS t)]))

(defn set-xyzt-units [^NiftiHeader h xyz t]
  (edit h #(set! (.xyztUnits ^NiftiHeader %)
                 (int (+ (kw->code Codes/UNITS (or xyz :unknown))
                         (kw->code Codes/UNITS (or t :unknown)))))))

(defn dim-info
  "Which voxel axes carry the frequency, phase and slice encoding directions,
  as a vector of three axis indices or nils."
  [^NiftiHeader h]
  (let [info (.dimInfo h)]
    (mapv #(let [v (bit-and (bit-shift-right info %) 3)]
             (when (pos? v) (dec v)))
          [0 2 4])))

(defn set-dim-info [^NiftiHeader h freq phase slice]
  (doseq [v [freq phase slice]]
    (when-not (or (nil? v) (<= 0 v 2))
      (throw (NiftiError. (str "dim_info axes must be 0, 1, 2 or nil, got " v)))))
  (edit h #(set! (.dimInfo ^NiftiHeader %)
                 (int (reduce bit-or 0
                              (map (fn [v shift] (if v (bit-shift-left (inc v) shift) 0))
                                   [freq phase slice] [0 2 4]))))))

(defn intent
  "Statistical intent as `[code-keyword name [p1 p2 p3]]`."
  [^NiftiHeader h]
  [(code->kw Codes/INTENT (.intentCode h))
   (.intentName h)
   [(.intentP1 h) (.intentP2 h) (.intentP3 h)]])

(defn slice-code [^NiftiHeader h] (code->kw Codes/SLICE (.sliceCode h)))

(defn descrip [^NiftiHeader h] (.descrip h))

(defn set-descrip [^NiftiHeader h s]
  (edit h #(set! (.descrip ^NiftiHeader %) (str s))))

;; ------------------------------------------------------------------ inspection

(defn ->map
  "Every header field as a map, for inspection and diffing."
  [^NiftiHeader h]
  {:version (.version h)
   :byte-order (if (big-endian? h) :big :little)
   :single-file? (.singleFile h)
   :magic (.magic h)
   :dim (vec (.dim h))
   :shape (shape h)
   :datatype (data-dtype h)
   :bitpix (.bitpix h)
   :pixdim (vec (.pixdim h))
   :zooms (zooms h)
   :vox-offset (.voxOffset h)
   :scl-slope (.sclSlope h)
   :scl-inter (.sclInter h)
   :qform-code (qform-code h)
   :sform-code (sform-code h)
   :quatern [(.quaternB h) (.quaternC h) (.quaternD h)]
   :qoffset [(.qoffsetX h) (.qoffsetY h) (.qoffsetZ h)]
   :srow (la/rows (.srow h))
   :xyzt-units (xyzt-units h)
   :dim-info (dim-info h)
   :intent (intent h)
   :slice-code (slice-code h)
   :slice-start (.sliceStart h)
   :slice-end (.sliceEnd h)
   :slice-duration (.sliceDuration h)
   :toffset (.toffset h)
   :cal-min (.calMin h)
   :cal-max (.calMax h)
   :descrip (.descrip h)
   :aux-file (.auxFile h)
   :extensions (mapv (fn [^Extension e] {:code (.code e) :size (alength (.content e))})
                     (.extensions h))})

(defn describe
  "A short human-readable summary of `h`."
  [^NiftiHeader h]
  (let [m (->map h)]
    (str/join "\n"
              [(format "NIfTI-%d %s %s" (:version m) (name (:byte-order m))
                       (if (:single-file? m) "single file" "hdr/img pair"))
               (format "  shape      %s" (:shape m))
               (format "  zooms      %s %s" (:zooms m) (name (first (:xyzt-units m))))
               (format "  datatype   %s (bitpix %d)" (name (:datatype m)) (:bitpix m))
               (format "  scaling    slope %s inter %s" (:scl-slope m) (:scl-inter m))
               (format "  qform      %s" (name (:qform-code m)))
               (format "  sform      %s" (name (:sform-code m)))
               "  affine"
               (str/join "\n" (map #(format "    %s" %) (best-affine h)))
               (format "  descrip    %s" (pr-str (:descrip m)))])))

(defn to-bytes
  "Serialise `h` plus its extension blocks."
  ^bytes [^NiftiHeader h]
  (.toBytesWithExtensions h))

(defn from-bytes
  "Parse a header (and any extensions) from the head of `buf`."
  ^NiftiHeader [^bytes buf]
  (NiftiHeader/read buf))
