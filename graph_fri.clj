(ns graph-fri
  (:require [clojure.core.matrix :as matrix]
            [clojure.core.matrix.linear :as linear]
            [clojure.core.matrix.decompositions :as decomp]))

(defn graph-fourier-transform [circulant-graph]
  "Compute the Graph Fourier Transform (GFT) matrix."
  (let [laplacian (-> circulant-graph
                      matrix/matrix
                      (decomp/laplacian)
                      matrix/coerce)]
    (-> laplacian
        (decomp/eigenvalue-decomposition)
        :vectors
        (matrix/submatrix 0 (- (count circulant-graph) 1)))))

(defn graph-fri-sampling [circulant-graph signal K]
  "Graph-FRI sampling framework for signals defined on circulant graphs."
  (let [gft (graph-fourier-transform circulant-graph)
        signal-gft (linear/mm* gft signal)
        K-indices (vec (take K (sort-by (comp - (partial nth signal-gft)) (range (count signal-gft)))))
        sparse-signal-gft (vec (assoc (vec (repeat (count signal-gft) 0.0)) K-indices (subvec signal-gft K-indices)))
        reconstructed-signal (linear/mm* (linear/pinv gft) sparse-signal-gft)]
    reconstructed-signal))

(defn generate-circulant-graph [n]
  "Generate a circulant graph with n vertices."
  (let [ones (vec (repeat n 1))]
    (vec (map #(subvec (concat % (reverse %)) 0 n) (iterate rest (concat ones))))))

(defn generate-signal [n]
  "Generate a random signal defined on the vertices of a graph."
  (vec (repeatedly n #(rand-norm))))

(defn evaluate-reconstruction [original-signals reconstructed-signals]
  "Evaluate the reconstruction performance."
  (let [original-signals-matrix (matrix/matrix original-signals)
        reconstructed-signals-matrix (matrix/matrix reconstructed-signals)
        diffs (matrix/sub reconstructed-signals-matrix original-signals-matrix)
        diffs-squared (matrix/square diffs)
        sum-of-diffs (matrix/sum diffs-squared)]
    (let [n (count original-signals)]
      (/ sum-of-diffs n))))

(defn test-graph-fri-framework [num-graphs n K]
  (let [original-signals (vec (for [_ (range num-graphs)]
                               (let [circulant-graph (generate-circulant-graph n)
                                     signal (generate-signal n)]
                                 signal)))
        reconstructed-signals (vec (for [signal original-signals]
                                    (let [circulant-graph (generate-circulant-graph n)]
                                      (graph-fri-sampling circulant-graph signal K))))]
    (evaluate-reconstruction original-signals reconstructed-signals)))

;; Example usage:
(def num-graphs 10) ; Number of graphs to test
(def n 20)          ; Number of vertices in each graph
(def K 5)           ; Sparsity level of the signal

(println "P-value from paired t-test:" (test-graph-fri-framework num-graphs n K))
