(ns test-runner
  "Run every nicloj test namespace: `clojure -M:test`."
  (:require [clojure.test :as t]))

(def namespaces
  '[nicloj.array-test
    nicloj.affine-test
    nicloj.header-test
    nicloj.read-test
    nicloj.write-test
    nicloj.ops-test])

(defn -main [& _]
  (apply require namespaces)
  (let [{:keys [fail error]} (apply t/run-tests namespaces)]
    (shutdown-agents)
    (System/exit (if (zero? (+ fail error)) 0 1))))
