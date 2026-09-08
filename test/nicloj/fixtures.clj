(ns nicloj.fixtures
  "Shared helpers for the test suite: the nibabel-generated corpus and
  approximate comparison."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]))

(def data-dir (io/file "test-data"))
(def out-dir (io/file data-dir "out"))

(def manifest
  "What nibabel reports for every file in the corpus, from
  `scripts/gen_testdata.py`."
  (delay
    (let [f (io/file data-dir "manifest.edn")]
      (when-not (.exists f)
        (throw (ex-info "test corpus missing; run scripts/gen_testdata.py first"
                        {:expected (str f)})))
      (edn/read-string (slurp f)))))

(defn corpus-file
  "Absolute path of a corpus file by name."
  [name]
  (.getPath (io/file data-dir name)))

(defn out-file
  "Path under test-data/out/ for files the tests write."
  [name]
  (.mkdirs out-dir)
  (.getPath (io/file out-dir name)))

(defn close?
  "Compare two numbers with a relative tolerance, treating NaN as equal."
  ([a b] (close? a b 1e-6))
  ([a b tol]
   (let [a (double a) b (double b)]
     (or (and (Double/isNaN a) (Double/isNaN b))
         (<= (Math/abs (- a b)) (* tol (max 1.0 (Math/abs a) (Math/abs b))))))))

(defn all-close?
  "Compare nested sequences of numbers elementwise."
  ([a b] (all-close? a b 1e-6))
  ([a b tol]
   (cond
     (and (number? a) (number? b)) (close? a b tol)
     (and (sequential? a) (sequential? b))
     (and (= (count a) (count b)) (every? true? (map #(all-close? %1 %2 tol) a b)))
     :else (= a b))))
