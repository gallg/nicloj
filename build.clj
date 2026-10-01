(ns build
  "Build tasks: `clojure -T:build javac`, `clojure -T:build clean`.

  Plain javac works just as well and needs no dependencies:
    javac -d target/classes $(find nicloj -name '*.java')"
  (:require [clojure.tools.build.api :as b]))

(def version
  "nicloj's version. Bump it here, and only here, with each release."
  "0.1.0")

(def class-dir "target/classes")

(defn javac
  "Compile every Java module under nicloj/ into target/classes."
  [_]
  (b/javac {:src-dirs ["nicloj"]
            :class-dir class-dir
            :javac-opts ["-Xlint:all" "-proc:none"]})
  (println "compiled Java modules ->" class-dir))

(defn clean [_]
  (b/delete {:path "target"})
  (println "removed target/"))
