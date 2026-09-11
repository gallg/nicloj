(ns nicloj.core.image
  "The in-memory NIfTI image: a header, a voxel-to-world affine and voxel data.

  Voxel data is held as it appears on disk and is loaded lazily; `fdata`
  applies the header's `scl_slope`/`scl_inter` to give real-world values, the
  way nibabel's `get_fdata` does."
  (:require [nicloj.core.header :as hdr]
            [nicloj.core.ndarray :as nd])
  (:import (nicloj.array NdArray)
           (nicloj.header NiftiHeader)))

(defrecord NiftiImage [^NiftiHeader header affine data])

(defn image?
  [x] (instance? NiftiImage x))

;; ------------------------------------------------------------------ accessors

(defn header ^NiftiHeader [img] (:header img))
(defn affine [img] (:affine img))
(defn shape [img] (hdr/shape (header img)))
(defn ndim [img] (count (shape img)))
(defn zooms [img] (hdr/zooms (header img)))
(defn data-dtype [img] (hdr/data-dtype (header img)))

(defn raw-data
  "Voxel values exactly as stored, before `scl_slope`/`scl_inter`, and in the
  file's own type, so an `int16` image costs two bytes a voxel."
  ^NdArray [img]
  (force (:data img)))

(defn fdata
  "Voxel values in real-world units, i.e. raw values with the header's scaling
  applied.

  Every element reads as a `double`, but the array is only *stored* as doubles
  when it has to be: without scaling this hands back `raw-data` untouched, in
  the file's own type. A scaled image widens, since the scaled values are no
  longer representable in it."
  ^NdArray [img]
  (let [[slope inter] (or (hdr/slope-inter (header img)) [1.0 0.0])]
    (nd/scaled (raw-data img) slope inter)))

(defn voxel
  "The real-world value at voxel index `idx`."
  ^double [img idx]
  (nd/value (fdata img) idx))

(defn loaded?
  "True when the voxel data is already in memory."
  [img]
  (let [d (:data img)]
    (or (not (delay? d)) (realized? d))))

;; ---------------------------------------------------------------- construction

(defn image
  "Build an image from `data` and a 4x4 `affine`.

  `data` may be an `NdArray` or row-major nested sequences. Options:

    :header  a template header to inherit non-geometric fields from
    :dtype   on-disk datatype; defaults to the header's, or `:float32`
    :version 1 (default) or 2, ignored when :header is given

  The affine is written into the header as nicloj writes new images: sform
  `:aligned` plus qform values coded `:unknown`. `data` is taken to be in
  real-world units, so header scaling is reset to the identity."
  [data affine & {:keys [header dtype version]}]
  (let [arr (nd/coerce data)
        base (or header (hdr/new-header :version (or version 1)))
        h (cond-> (hdr/clear-scaling (hdr/set-shape base (nd/shape arr)))
            dtype (hdr/set-data-dtype dtype)
            affine (hdr/set-affine affine))]
    (->NiftiImage h (or affine (hdr/best-affine h)) arr)))

(defn from-parts
  "Assemble an image from an already-prepared header and voxel data. `data` may
  be a delay, so that reading voxels can be postponed."
  [^NiftiHeader h data]
  (->NiftiImage h (hdr/best-affine h) data))

(defn with-data
  "Replace the voxel data, updating `dim` to match. The new values are taken to
  be in real-world units, so header scaling is reset."
  [img data]
  (let [arr (nd/coerce data)]
    (-> img
        (assoc :header (hdr/clear-scaling (hdr/set-shape (header img) (nd/shape arr))))
        (assoc :data arr))))

(defn with-affine
  "Replace the affine, writing it into the header's sform and qform."
  [img affine]
  (-> img
      (assoc :header (hdr/set-affine (header img) affine))
      (assoc :affine affine)))

(defn with-header
  "Replace the header, keeping the current voxel data and the values it stands
  for. The affine is re-read from the new header, but `scl_slope`/`scl_inter`
  carry over from the old one, since the stored data is relative to those."
  [img ^NiftiHeader h]
  (let [[slope inter] (or (hdr/slope-inter (header img)) [1.0 0.0])
        h (-> h (hdr/set-shape (shape img)) (hdr/set-slope-inter slope inter))]
    (->NiftiImage h (hdr/best-affine h) (:data img))))

(defn with-dtype
  "Change the on-disk datatype used when the image is next written."
  [img dtype]
  (assoc img :header (hdr/set-data-dtype (header img) dtype)))

;; ---------------------------------------------------------------------- print

(defn- summary [img]
  (format "#nicloj/image{:shape %s :dtype %s :zooms %s :loaded? %s}"
          (shape img) (data-dtype img) (mapv float (zooms img)) (loaded? img)))

(defmethod print-method NiftiImage [img ^java.io.Writer w]
  (.write w ^String (summary img)))

(defn describe
  "A multi-line summary of the image and its header."
  [img]
  (str (summary img) "\n" (hdr/describe (header img))))
