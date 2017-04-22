;;;; package.lisp

(defpackage #:noise-solver
  (:use #:cl #:alexandria)
  (:export :collect-numeric-directories
           :average-data
           :average-momentum
           :generate-noise
           :noise-solver
           :main))
