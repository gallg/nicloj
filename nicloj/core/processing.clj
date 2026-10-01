(ns nicloj.core.processing
  "Resampling images onto other voxel grids, after `nibabel.processing`.

  Interpolation is nearest-neighbour (`:order 0`) or linear (`:order 1`, the
  default), with scipy's `mode='constant'` edges: a point outside the source
  grid takes `:cval`. nibabel defaults to cubic splines, which nicloj does not
  have yet."
  (:require [nicloj.core.image :as img]
            [nicloj.core.linalg :as la])
  (:import (nicloj.array Resample)
           (nicloj.header NiftiError)))

(defn- adapt-affine
  "Embed a 4x4 affine in the `(n+1)`-square one for `n` voxel axes, the axes
  past the third mapping to themselves, as nibabel's `adapt_affine` does."
  [affine n]
  (vec (for [r (range (inc n))]
         (vec (for [c (range (inc n))]
                (cond
                  (and (< r 3) (< c 3)) (double (get-in affine [r c]))
                  (and (< r 3) (= c n)) (double (get-in affine [r 3]))
                  (= r c) 1.0
                  :else 0.0))))))

(defn resample-from-to
  "Resample `from` onto the voxel grid of `to`, which is an image or a
  `[shape affine]` pair.

  `from` must be at least 3D. Axes past the third are carried through
  unchanged; when `to` has fewer axes than `from`, the missing ones are taken
  from `from`, so a 4D series can be resampled onto a 3D template. Options:

    :order  0 nearest neighbour, 1 linear (default)
    :cval   value for points outside `from`, default 0.0

  The result keeps `from`'s header fields and voxel type, as nibabel's does."
  [from to & {:keys [order cval] :or {order 1 cval 0.0}}]
  (let [[to-shape to-affine] (if (img/image? to) [(img/shape to) (img/affine to)] to)
        from-shape (img/shape from)
        n (count from-shape)
        _ (when (< n 3)
            (throw (NiftiError. (str "from image must be at least 3D, got shape " from-shape))))
        _ (when (or (< (count to-shape) 3) (> (count to-shape) n))
            (throw (NiftiError. (str "cannot resample shape " from-shape " onto shape " to-shape))))
        to-shape (into (vec to-shape) (drop (count to-shape) from-shape))
        to-vox->from-vox (la/mmul (la/minverse (adapt-affine (img/affine from) n))
                                  (adapt-affine to-affine n))
        data (Resample/affine (img/fdata from) (la/arr to-vox->from-vox)
                              (int-array to-shape) (int order) (double cval))]
    (img/image data to-affine :header (img/header from))))

(defn resample-to-output
  "Resample a 3D image onto an axis-aligned RAS grid that just covers it, with
  `:voxel-sizes` (a number or three; default 1) as the new voxel edges. Takes
  the same `:order` and `:cval` as `resample-from-to`."
  [image & {:keys [voxel-sizes order cval] :or {voxel-sizes 1.0 order 1 cval 0.0}}]
  (let [shape (img/shape image)
        _ (when-not (= 3 (count shape))
            (throw (NiftiError. (str "resample-to-output needs a 3D image, got shape " shape))))
        sizes (mapv double (if (number? voxel-sizes) (repeat 3 voxel-sizes) voxel-sizes))
        _ (when-not (and (= 3 (count sizes)) (every? pos? sizes))
            (throw (NiftiError. (str "need three positive voxel sizes, got " voxel-sizes))))
        [x y z] (mapv dec shape)
        corners (for [i [0 x] j [0 y] k [0 z]]
                  (la/apply-affine (img/affine image) [i j k]))
        lo (apply mapv min corners)
        hi (apply mapv max corners)
        out-shape (mapv #(long (inc (Math/ceil (/ (- %2 %1) %3)))) lo hi sizes)
        out-affine (-> (la/identity 4)
                       (assoc-in [0 0] (sizes 0)) (assoc-in [1 1] (sizes 1)) (assoc-in [2 2] (sizes 2))
                       (assoc-in [0 3] (lo 0)) (assoc-in [1 3] (lo 1)) (assoc-in [2 3] (lo 2)))]
    (resample-from-to image [out-shape out-affine] :order order :cval cval)))
