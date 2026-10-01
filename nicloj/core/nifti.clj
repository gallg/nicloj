(ns nicloj.core.nifti
  "Reading and writing NIfTI-1 and NIfTI-2 files.

  Both layouts are handled in either byte order, as a single `.nii` file or a
  `.hdr`/`.img` pair, plain or gzipped. Header extensions are carried through
  unchanged."
  (:refer-clojure :exclude [load])
  (:require [nicloj.core.header :as hdr]
            [nicloj.core.image :as img]
            [nicloj.core.linalg :as la]
            [nicloj.core.ndarray :as nd]
            [nicloj.io.paths :as paths]
            [nicloj.io.stream :as stream])
  (:import (nicloj.array Codec NdArray)
           (nicloj.header DataType NiftiError NiftiHeader)))

(defn- concat-bytes ^bytes [^bytes a ^bytes b]
  (let [out (byte-array (+ (alength a) (alength b)))]
    (System/arraycopy a 0 out 0 (alength a))
    (System/arraycopy b 0 out (alength a) (alength b))
    out))

(defn- checked-count ^long [^long n what]
  (when (> n Integer/MAX_VALUE)
    (throw (NiftiError. (str what " is " n " bytes, over the 2 GiB nicloj can hold in one array"))))
  n)

;; --------------------------------------------------------------------- reading

(defn- read-struct
  "Parse the header of the NIfTI file at `path`.

  Returns `[header data-offset]`, where the offset is into `path` for a single
  file and into the companion `.img` for a pair."
  [path single?]
  (with-open [in (stream/input-stream path)]
    (let [head (stream/read-n! in 4)
          version (aget (NiftiHeader/probe head) 0)
          size (if (= 2 version) NiftiHeader/NIFTI2_SIZE NiftiHeader/NIFTI1_SIZE)
          fixed (concat-bytes head (stream/read-n! in (- size 4)))
          probed (NiftiHeader/read fixed)
          offset (if single? (max (.voxOffset probed) (+ size 4)) (.voxOffset probed))
          _ (when (neg? offset) (throw (NiftiError. (str "negative vox_offset " offset))))
          extra (if single?
                  (.readNBytes in (int (checked-count (- offset size) "header extensions")))
                  (.readAllBytes in))]
      [(NiftiHeader/read (concat-bytes fixed extra)) offset])))

(defn- read-voxels
  "Read the voxel block of `h` from `path`, starting at `offset`. Values come
  back exactly as stored, without `scl_slope`/`scl_inter`."
  ^NdArray [path offset ^NiftiHeader h]
  (with-open [in (stream/input-stream path)]
    (stream/skip-n! in offset)
    (Codec/decode (stream/read-n! in (checked-count (.dataBytes h) "voxel data"))
                  (.shape h) (.dataType h) (.order h) 1.0 0.0)))

(defn read-header
  "Read just the header of the NIfTI image named by `path`."
  ^NiftiHeader [path]
  (let [{:keys [header single?]} (paths/existing-image path)]
    (first (read-struct header single?))))

(defn load
  "Read the NIfTI image named by `path`.

  `path` may name a `.nii[.gz]` file or either half of a `.hdr`/`.img[.gz]`
  pair. Voxel data is read on first use; pass `:eager? true` to read it now."
  [path & {:keys [eager?]}]
  (let [{hdr-path :header img-path :image :keys [single?]} (paths/existing-image path)
        [h offset] (read-struct hdr-path single?)
        data (delay (read-voxels img-path offset h))]
    (img/from-parts h (if eager? (force data) data))))

(defn load-slab
  "Read only indices `start` (inclusive) to `stop` (exclusive) of the last
  axis of the image named by `path` -- volumes of a 4D image, slices of a 3D
  one. Just those bytes are read: an uncompressed file is skipped past, a
  gzipped one streamed past.

  The header keeps its scaling and xform codes, with the forms shifted so
  voxels keep their world positions. An image with neither form set gets the
  shifted fallback affine as an `:aligned` sform when the fallback for its new
  shape would put them elsewhere, as nibabel's slicer does."
  [path start stop]
  (let [{hdr-path :header img-path :image :keys [single?]} (paths/existing-image path)
        [^NiftiHeader h offset] (read-struct hdr-path single?)
        shape (hdr/shape h)
        axis (dec (count shape))]
    (when-not (and (integer? start) (integer? stop) (<= 0 start stop (nth shape axis)))
      (throw (NiftiError. (str "slab [" start " " stop ") out of range for shape " shape))))
    (let [per (.dataBytes ^NiftiHeader (hdr/set-shape h (assoc shape axis 1)))
          shift (cond-> (la/identity 4) (< axis 3) (assoc-in [axis 3] (double start)))
          sub (hdr/set-shape h (assoc shape axis (- stop start)))
          moved (la/mmul (hdr/best-affine h) shift)
          sub (cond-> sub
                (and (< axis 3) (not= 0 (.sformCode h)))
                (hdr/set-sform (la/mmul (hdr/sform h) shift))

                (and (< axis 3) (not= 0 (.qformCode h)))
                (hdr/set-qform (la/mmul (hdr/qform h) shift))

                (and (= 0 (.sformCode h) (.qformCode h))
                     (not (la/close? moved (hdr/base-affine sub) 1e-9)))
                (hdr/set-sform moved :aligned))]
      (img/from-parts sub (read-voxels img-path (+ offset (* start per)) sub)))))

(defn from-bytes
  "Read a complete single-file NIfTI image from a byte array."
  [^bytes buf]
  (let [h (NiftiHeader/read buf)
        offset (max (.voxOffset h) (+ (.structSize h) 4))
        end (+ offset (.dataBytes h))]
    ;; copyOfRange pads with zeros past the end, so check before slicing rather
    ;; than handing back a silently truncated volume.
    (when (> end (alength buf))
      (throw (NiftiError. (str "NIfTI byte array is truncated: voxel data ends at "
                               end " but the array is " (alength buf) " bytes"))))
    (img/from-parts h (Codec/decode (java.util.Arrays/copyOfRange buf (int offset) (int end))
                                    (.shape h) (.dataType h) (.order h) 1.0 0.0))))

;; --------------------------------------------------------------------- writing

(defn- target-dtype ^DataType [img dtype]
  (DataType/fromLabel (name (hdr/dtype (or dtype (hdr/data-dtype (img/header img)))))))

(defn- resolve-scaling
  "Decide what to write: the voxel array, the `scl_slope`/`scl_inter` to record
  in the header, and whether the array is already in stored (raw) units."
  [img ^DataType dt scaling version]
  (let [h (img/header img)]
    (case scaling
      :keep (do
              (when-not (= (.datatype h) (.code dt) (.code (.dtype (img/raw-data img))))
                (throw (NiftiError. ":scaling :keep cannot change the datatype")))
              (let [[slope inter] (or (hdr/slope-inter h) [1.0 0.0])]
                {:array (img/raw-data img) :slope slope :inter inter :raw? true}))
      :auto (let [arr (img/fdata img)
                  [slope inter] (vec (Codec/autoScale arr dt (= 1 version)))]
              {:array arr :slope slope :inter inter})
      (let [[slope inter] (when (and (sequential? scaling) (= 2 (count scaling))) scaling)]
        (when-not (and (number? slope) (number? inter))
          (throw (NiftiError. (str ":scaling must be :auto, :keep or [slope inter], got " scaling))))
        {:array (img/fdata img) :slope (double slope) :inter (double inter)}))))

(defn- prepare
  "A copy of the image's header, updated to describe what is about to be written."
  ^NiftiHeader [img ^NdArray arr ^DataType dt slope inter version single?]
  (let [h (.copy (img/header img))]
    (set! (.version h) (int version))
    (set! (.singleFile h) (boolean single?))
    (.setShape h (int-array (nd/shape arr)))
    (.setDataType h dt)
    (set! (.sclSlope h) (double slope))
    (set! (.sclInter h) (double inter))
    (set! (.voxOffset h) (long (if single? (.minVoxOffset h) 0)))
    h))

(defn- encode
  "Serialise `img` into `[header-bytes voxel-bytes]`."
  [img single? {:keys [dtype scaling version] :or {scaling :auto}}]
  (let [version (or version (.version (img/header img)))
        _ (when-not (#{1 2} version)
            (throw (NiftiError. (str ":version must be 1 or 2, got " version))))
        dt (target-dtype img dtype)
        {:keys [array slope inter raw?]} (resolve-scaling img dt scaling version)
        [slope inter] (if (= 1 version)
                        [(double (unchecked-float slope)) (double (unchecked-float inter))]
                        [slope inter])
        _ (when-not (and (Double/isFinite slope) (not (zero? slope)) (Double/isFinite inter))
            (throw (NiftiError. (str "cannot write scl_slope " slope ", scl_inter " inter
                                     ": both must be finite and the slope non-zero"))))
        h (prepare img array dt slope inter version single?)]
    [(.toBytesWithExtensions h)
     (if raw?
       (Codec/encode array dt (.order h) 1.0 0.0)
       (Codec/encode array dt (.order h) slope inter))]))

(defn save
  "Write `img` to `path` and return `path`.

  The name decides the layout: `.nii[.gz]` for one file, `.hdr`/`.img[.gz]`
  for a pair, a `.gz` suffix for gzip. Options:

    :dtype    on-disk datatype, defaulting to the image header's
    :version  1 or 2, defaulting to the image header's
    :scaling  `:auto` (default) chooses `scl_slope`/`scl_inter` so the values
              survive the target datatype; `:keep` reuses the header's scaling
              and writes the raw voxels back untouched, which round-trips a
              loaded file byte for byte; a `[slope inter]` pair sets them"
  [img path & {:as opts}]
  (let [{hdr-path :header img-path :image :keys [single?]} (paths/resolve-image path)
        [head body] (encode img single? opts)]
    (if single?
      (stream/write-all! hdr-path head body)
      (do (stream/write-all! hdr-path head)
          (stream/write-all! img-path body)))
    path))

(defn ->bytes
  "Serialise `img` as a complete single-file NIfTI image."
  ^bytes [img & {:as opts}]
  (apply concat-bytes (encode img true opts)))
