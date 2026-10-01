(ns nicloj.io.paths
  "Resolving NIfTI file names.

  A NIfTI image lives either in one `.nii` file or in a `.hdr`/`.img` pair,
  each optionally gzipped. Given any one of those names, `resolve-image`
  returns every path involved."
  (:require [clojure.string :as str])
  (:import (java.io File)))

(defn gz?
  "True when `path` carries a `.gz` suffix."
  [path]
  (str/ends-with? (str/lower-case (str path)) ".gz"))

(defn- split-gz
  "Split `path` into its non-gzip stem and a `.gz` suffix (possibly empty)."
  [path]
  (if (gz? path)
    [(subs path 0 (- (count path) 3)) (subs path (- (count path) 3))]
    [path ""]))

(defn- cut-ext
  "Split `s` into everything but its last four characters and those characters,
  lower-cased. A shorter string yields an empty stem, which no extension matches."
  [s]
  (let [at (max 0 (- (count s) 4))]
    [(subs s 0 at) (str/lower-case (subs s at))]))

(defn resolve-image
  "Locate the files making up the NIfTI image named by `path`.

  Returns `{:header p, :image p, :single? bool}`. For a `.nii` image both
  paths are the same file; for a pair, `:image` mirrors the gzip suffix of the
  name given, and the companion's extension is upper case when the given one
  is, as in nibabel. Throws when the extension is not a NIfTI one."
  [path]
  (let [path (str path)
        [stem gz] (split-gz path)
        [base ext] (cut-ext stem)
        given (subs stem (count base))
        companion (fn [e] (str base (if (= given (str/upper-case given)) (str/upper-case e) e) gz))]
    (case ext
      ".nii" {:header path :image path :single? true}
      ".hdr" {:header path :image (companion ".img") :single? false}
      ".img" {:header (companion ".hdr") :image path :single? false}
      (throw (ex-info (str "not a NIfTI filename: " path
                           " (expected .nii, .nii.gz, .hdr, .img or a .gz of those)")
                      {:path path})))))

(defn- on-disk
  "`path` itself, or its gzip sibling when only that one exists."
  [path]
  (if (.exists (File. ^String path))
    path
    (let [[stem gz] (split-gz path)
          alt (if (seq gz) stem (str path ".gz"))]
      (if (.exists (File. ^String alt)) alt path))))

(defn existing-image
  "Like `resolve-image`, but tolerant of a pair whose two halves disagree about
  gzip: each half falls back to the sibling that is actually on disk."
  [path]
  (let [{:keys [single?] :as m} (resolve-image path)]
    (if single?
      m
      (-> m (update :header on-disk) (update :image on-disk)))))
