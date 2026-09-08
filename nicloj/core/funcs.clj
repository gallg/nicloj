(ns nicloj.core.funcs
  "Everyday operations on whole images: cropping, splitting, joining and
  mapping between voxel and world coordinates."
  (:require [nicloj.core.image :as img]
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

  `specs` holds one entry per axis: nil for the whole axis, or
  `[start stop]` / `[start stop step]` for a range. Steps must be positive.
  Only the first three axes affect the affine."
  [image specs]
  (let [shape (img/shape image)
        _ (when (> (count specs) (count shape))
            (throw (NiftiError. (str "got " (count specs) " slice specs for shape " shape))))
        specs (mapv (fn [axis spec]
                      (let [[start stop step] spec]
                        [(or start 0) (or stop (nth shape axis)) (or step 1)]))
                    (range (count shape))
                    (concat specs (repeat nil)))
        transform (reduce (fn [m axis]
                            (let [[start _ step] (nth specs axis)]
                              (-> m
                                  (assoc-in [axis axis] (double step))
                                  (assoc-in [axis 3] (double start)))))
                          (la/identity 4)
                          (range (min 3 (count shape))))]
    (img/image (nd/slice (img/fdata image) specs)
               (la/mmul (img/affine image) transform)
               :header (img/header image))))
