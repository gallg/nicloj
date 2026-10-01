(ns nicloj.core.orientation
  "Axis orientations: which way round an image's voxel axes run, and how to
  turn one arrangement into another.

  An orientation is a vector with one `[out-axis flip]` entry per voxel axis;
  `flip` is `1` or `-1`, and an entry of `[nil nil]` marks a voxel axis with no
  matching world axis. Axis codes are the usual `\"R\" \"A\" \"S\"` letters."
  (:require [nicloj.core.header :as hdr]
            [nicloj.core.image :as img]
            [nicloj.core.linalg :as la]
            [nicloj.core.ndarray :as nd])
  (:import (nicloj.affine Orientations)
           (nicloj.header NiftiError)))

(def canonical
  "The orientation of RAS-ordered data: no transpose, no flips."
  [[0 1] [1 1] [2 1]])

(defn- ->arr ^"[[D" [ornt]
  (la/arr (map (fn [[axis flip]]
                 [(if (nil? axis) Double/NaN (double axis))
                  (if (nil? flip) Double/NaN (double flip))])
               ornt)))

(defn- ->ornt [^"[[D" a]
  (mapv (fn [^doubles row]
          [(when-not (Double/isNaN (aget row 0)) (long (aget row 0)))
           (when-not (Double/isNaN (aget row 1)) (long (aget row 1)))])
        a))

(defn io-orientation
  "How each input axis of `affine` maps onto its output axes. Pass `tol` to
  override the singular-value cutoff used to detect dropped axes."
  ([affine] (io-orientation affine 0))
  ([affine tol] (->ornt (Orientations/ioOrientation (la/arr affine) (double tol)))))

(defn ornt->axcodes
  "Axis labels for an orientation, with nil for dropped axes."
  [ornt]
  (vec (Orientations/toAxcodes (->arr ornt) Orientations/RAS_LABELS)))

(defn axcodes->ornt
  "The orientation described by axis labels such as `[\"R\" \"A\" \"S\"]`."
  [codes]
  (->ornt (Orientations/fromAxcodes (into-array String (map #(when % (name %)) codes))
                                    Orientations/RAS_LABELS)))

(defn axcodes
  "Axis labels of an affine's voxel axes, e.g. `[\"L\" \"A\" \"S\"]`."
  [affine]
  (ornt->axcodes (io-orientation affine)))

(defn ornt-transform
  "The orientation that takes data from `start` to `end`."
  [start end]
  (->ornt (Orientations/transform (->arr start) (->arr end))))

(defn inv-ornt-aff
  "The affine undoing the flips and transpose of `ornt` on data of `shape`, so
  that `(mmul old-affine (inv-ornt-aff ornt shape))` is the new affine."
  [ornt shape]
  (la/rows (Orientations/invOrntAff (->arr ornt) (int-array shape))))

(defn apply-orientation
  "Flip and transpose a voxel array according to `ornt`."
  [arr ornt]
  (let [a (nd/coerce arr)
        _ (when (> (count ornt) (nd/ndim a))
            (throw (NiftiError. (str "orientation has " (count ornt) " axes but the data only "
                                     (nd/ndim a)))))
        flipped (reduce (fn [acc [axis [_ flip]]]
                          (if (== -1 flip) (nd/flip acc axis) acc))
                        a
                        (map-indexed vector ornt))
        perm (vec (Orientations/axisPermutation (->arr ornt)))]
    (nd/transpose flipped (into perm (range (count ornt) (nd/ndim flipped))))))

(defn as-reoriented
  "Rearrange `img`'s voxel axes according to `ornt`, adjusting the affine and
  `dim_info` so the image still describes the same physical volume."
  [img ornt]
  (if (= (vec ornt) canonical)
    img
    (let [data (apply-orientation (img/fdata img) ornt)
          affine (la/mmul (img/affine img) (inv-ornt-aff ornt (img/shape img)))
          remap (fn [axis] (when axis (long (first (nth ornt axis)))))
          [freq phase slice] (hdr/dim-info (img/header img))
          h (hdr/set-dim-info (img/header img) (remap freq) (remap phase) (remap slice))]
      (img/image data affine :header h))))

(defn as-closest-canonical
  "Reorient `img` to the RAS arrangement closest to its current affine."
  [img]
  (as-reoriented img (io-orientation (img/affine img))))
