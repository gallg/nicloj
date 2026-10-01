(ns nicloj.core.funcs
  "Everyday operations on whole images: cropping, splitting, joining and
  mapping between voxel and world coordinates."
  (:require [nicloj.core.header :as hdr]
            [nicloj.core.image :as img]
            [nicloj.core.linalg :as la]
            [nicloj.core.ndarray :as nd])
  (:import (nicloj.header NiftiError)))

(defn voxel->world
  "World coordinates of a voxel index."
  [image idx]
  (la/apply-affine (img/affine image) idx))

(defn world->voxel
  "Voxel index of a world coordinate, as fractional numbers."
  [image xyz]
  (la/apply-affine (la/minverse (img/affine image)) xyz))

(defn squeeze-image
  "Drop trailing length-1 dimensions past the third, leaving the first three
  alone. Matches nibabel's `squeeze_image`."
  [image]
  (let [shape (img/shape image)]
    (if (or (< (count shape) 4) (not-any? #(= 1 %) (drop 3 shape)))
      image
      (let [trimmed (loop [s shape]
                      (if (and (> (count s) 3) (= 1 (peek s)))
                        (recur (pop s))
                        s))]
        (img/with-data image (nd/reshape (img/fdata image) trimmed))))))

(defn concat-images
  "Join images into one, taking the header and affine from the first.

  Without `:axis` the images are stacked on a new trailing axis; with an axis
  they are joined along that existing one. Affines must match unless
  `:check-affines? false`."
  [images & {:keys [axis check-affines?] :or {check-affines? true}}]
  (when (empty? images)
    (throw (NiftiError. "nothing to concatenate")))
  (let [affine (img/affine (first images))
        arrays (mapv img/fdata images)
        shape (nd/shape (first arrays))]
    (when check-affines?
      (doseq [i images]
        (when-not (la/close? (img/affine i) affine 0.0)
          (throw (NiftiError. "affines do not match")))))
    (when-not axis
      (doseq [[i a] (map-indexed vector arrays)]
        (when-not (= shape (nd/shape a))
          (throw (NiftiError. (str "image " i " has shape " (nd/shape a)
                                   ", not the first image's " shape))))))
    (let [parts (if axis arrays (mapv #(nd/reshape % (conj shape 1)) arrays))]
      (img/image (nd/concat parts (or axis (count shape))) affine
                 :header (img/header (first images))))))

(defn four-to-three
  "Split a 4D image into a vector of 3D images, one per volume."
  [image]
  (let [shape (img/shape image)]
    (when-not (= 4 (count shape))
      (throw (NiftiError. (str "expected a 4D image, got shape " shape))))
    (let [data (img/fdata image)
          volume (subvec shape 0 3)]
      (mapv (fn [i]
              (img/image (nd/reshape (nd/slice data [nil nil nil [i (inc i) 1]]) volume)
                         (img/affine image)
                         :header (img/header image)))
            (range (nth shape 3))))))

(defn slice-image
  "Crop or subsample `image`, adjusting the affine so the voxels keep their
  world positions.

  `specs` holds one entry per axis: nil for the whole axis,
  `[start stop]` / `[start stop step]` for a range, or, past the third axis,
  an integer, which picks one index and drops the axis, as nibabel's
  `img.slicer[..., 3]` does; negative integers count from the end. Steps must
  be positive. Only the first three axes affect the affine."
  [image specs]
  (let [shape (img/shape image)
        _ (when (> (count specs) (count shape))
            (throw (NiftiError. (str "got " (count specs) " slice specs for shape " shape))))
        given (vec (take (count shape) (concat specs (repeat nil))))
        kept (vec (remove #(integer? (given %)) (range (count shape))))
        specs (mapv (fn [axis spec]
                      (if (integer? spec)
                        (let [n (nth shape axis)
                              i (if (neg? spec) (+ n spec) spec)]
                          (when (< axis 3)
                            (throw (NiftiError. (str "an integer index would drop spatial axis " axis
                                                     "; use a one-voxel range [i (inc i)] instead"))))
                          (when-not (< -1 i n)
                            (throw (NiftiError. (str "index " spec " out of range for axis " axis
                                                     " of size " n))))
                          [i (inc i) 1])
                        (let [[start stop step] spec]
                          [(or start 0) (or stop (nth shape axis)) (or step 1)])))
                    (range (count shape))
                    given)
        transform (reduce (fn [m axis]
                            (let [[start _ step] (nth specs axis)]
                              (-> m
                                  (assoc-in [axis axis] (double step))
                                  (assoc-in [axis 3] (double start)))))
                          (la/identity 4)
                          (range (min 3 (count shape))))
        data (nd/slice (img/fdata image) specs)
        data (nd/reshape data (mapv (nd/shape data) kept))
        ;; Dropped axes take their voxel sizes with them.
        h (-> (img/header image)
              (hdr/set-shape (nd/shape data))
              (hdr/set-zooms (mapv (img/zooms image) kept)))]
    (img/image data (la/mmul (img/affine image) transform) :header h)))
