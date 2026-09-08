(ns nicloj.io.stream
  "Byte-level file access with transparent gzip.

  Reading sniffs the gzip magic number rather than trusting the file name, so
  a mislabelled `.nii` holding compressed bytes still loads. Writing
  compresses when the name ends in `.gz`."
  (:require [nicloj.io.paths :as paths])
  (:import (java.io BufferedInputStream BufferedOutputStream File FileInputStream
                    FileOutputStream InputStream OutputStream PushbackInputStream)
           (java.util.zip GZIPInputStream GZIPOutputStream)))

(def ^:private buffer-size 65536)

(defn- gzip-stream?
  "Peek at the first two bytes of `in` for the gzip magic 1f 8b."
  [^PushbackInputStream in]
  (let [b0 (.read in)]
    (if (neg? b0)
      false
      (let [b1 (.read in)]
        (when-not (neg? b1) (.unread in b1))
        (.unread in b0)
        (and (= 0x1f b0) (= 0x8b b1))))))

(defn input-stream
  "Open `path` for reading, decompressing if its contents are gzipped."
  ^InputStream [path]
  (let [raw (PushbackInputStream. (FileInputStream. (File. (str path))) 2)]
    (BufferedInputStream.
     (if (gzip-stream? raw) (GZIPInputStream. raw buffer-size) raw)
     buffer-size)))

(defn output-stream
  "Open `path` for writing, compressing when the name ends in `.gz`."
  ^OutputStream [path]
  (let [raw (FileOutputStream. (File. (str path)))]
    (if (paths/gz? path)
      (GZIPOutputStream. raw (int buffer-size))
      (BufferedOutputStream. raw buffer-size))))

(defn read-n!
  "Read exactly `n` bytes from `in`, or throw on a short read."
  ^bytes [^InputStream in n]
  (let [buf (byte-array n)]
    (loop [off 0]
      (when (< off n)
        (let [got (.read in buf off (- n off))]
          (if (neg? got)
            (throw (ex-info (str "unexpected end of file after " off " of " n " bytes")
                            {:read off :expected n}))
            (recur (+ off got))))))
    buf))

(defn skip-n!
  "Skip exactly `n` bytes of `in`, or throw if the stream ends first."
  [^InputStream in n]
  (loop [left (long n)]
    (when (pos? left)
      (let [got (.skip in left)]
        (if (pos? got)
          (recur (- left got))
          (if (neg? (.read in))
            (throw (ex-info (str "unexpected end of file while skipping " n " bytes")
                            {:bytes n}))
            (recur (dec left))))))))

(defn write-all!
  "Write the given byte arrays to `path` in order."
  [path & arrays]
  (with-open [out (output-stream path)]
    (doseq [^bytes a arrays] (.write out a))))
