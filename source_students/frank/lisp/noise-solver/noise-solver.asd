;;;; noise-solver.asd

(asdf:defsystem #:noise-solver
  :description "Scripts for noise_solver"
  :author "Frank Stuckenberg <frank.stuckenberg@uni-bremen.de>"
  :license "GPLv3"
  :serial t
  :depends-on (:uiop :ieee-floats :alexandria)
  :components ((:file "package")
               (:file "binary-data")
               (:file "noise-solver" :depends-on ("binary-data"))))
