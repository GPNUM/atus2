;;;; noise-solver.lisp

(in-package #:noise-solver)

(setf *read-default-float-format* 'double-float)

(defun time-slice (start end i n)
  (assert (>= end start 0))
  (assert (>= n i 0))
  (loop for j below i with x = start
     do (setf x (sqrt (+ (/ (- (expt end 2) (expt start 2)) n) (expt x 2))))
     finally (return x)))

(defun gen-gen.sh (start N)
		(with-open-file (out (format nil "gen~d.sh" N) :direction :output :if-exists :supersede :if-does-not-exist :create)
		  (format out "#!/bin/sh~%~%")
		  (format out "~&../run-noise-solver.lisp generate ~d ~d generate_noise.xml & ~%"  (+ start (* 10 0)) (+ 9 start (* 10 0)))
		  (format out "~&../run-noise-solver.lisp generate ~d ~d generate_noise.xml & ~%"  (+ start (* 10 1)) (+ 9 start (* 10 1)))
		  (format out "~&../run-noise-solver.lisp generate ~d ~d generate_noise.xml & ~%"  (+ start (* 10 2)) (+ 9 start (* 10 2)))
		  (format out "~&../run-noise-solver.lisp generate ~d ~d generate_noise.xml & ~%"  (+ start (* 10 3)) (+ 9 start (* 10 3)))
		  (format out "~&../run-noise-solver.lisp generate ~d ~d generate_noise.xml & ~%"  (+ start (* 10 4)) (+ 9 start (* 10 4)))
		  (format out "~&../run-noise-solver.lisp generate ~d ~d generate_noise.xml & ~%"  (+ start (* 10 5)) (+ 9 start (* 10 5)))
		  (format out "~&~%wait~%")))

(defun gen-all-gen.sh (start)
  (gen-gen.sh (+ start (* 60 0)) 1)
  (gen-gen.sh (+ start (* 60 1)) 2)
  (gen-gen.sh (+ start (* 60 2)) 3)
  (gen-gen.sh (+ start (* 60 3)) 4)
  (gen-gen.sh (+ start (* 60 4)) 5))

(defun gen-freeprop (filename start end)
  (with-open-file (out filename :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format out "#!/bin/bash
#PBS -j oe
#PBS -l walltime=11:59:00
#PBS -l nodes=~d:ppn=24
#PBS -l feature=mpp
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -A hbp00035

source $HOME/frank.bashrc
export MY_NO_OF_THREADS=8
export OMP_NUM_THREADS=8

cd wide~%~%" (/ (- end start -1) 3))
    (loop :for i :from start :to end :by 3
       :do (format out "~&aprun -n 1 -d 24 ../freeprop-helper.sh ~d ~d ~d &" i (+ i 1) (+ i 2)))
    (format out "~%~%wait~&")))

(defun generate-run.sh (noise mstart mend step g1 &optional (files 3))
  (loop :for file :from 1 :to files
     :for start = (* (round (time-slice mstart mend (1- file) files) step) step)
     :for end = (* (round (time-slice mstart mend file files) step) step)
     :do (format t "~a ~a ~a ~a~%" start end file files)
       (with-open-file (out (format nil "~a-~a-~a-g1-~a.sh" file 3 noise g1) :direction :output :if-exists :supersede :if-does-not-exist :create)
         (format out "#!/bin/sh~%~%")
         (format out "~&../run-noise-solver.lisp solve $1 $2 ~{if-~a-chirp-~a-g1-~a.xml ~}~%~%"
                 (loop
                    for i
                    from (round (time-slice start end 0 3) step)
                    below (round (time-slice start end 1 3) step)
                    append (list noise (+ 0 (* i step)) g1)))
         (format out "~&../run-noise-solver.lisp solve $1 $2 ~{if-~a-chirp-~a-g1-~a.xml ~}~%~%"
                 (loop
                    for i
                    from (round (time-slice start end 1 3) step)
                    below (round (time-slice start end 2 3) step)
                    append (list noise (+ 0 (* i step)) g1)))
         (format out "~&../run-noise-solver.lisp solve $1 $2 ~{if-~a-chirp-~a-g1-~a.xml ~}~%~%"
                 (loop
                    for i
                    from (round (time-slice start end 2 3) step)
                    below (if (= file files)
                              (1+ (round (time-slice start end 3 3) step))
                              (round (time-slice start end 3 3) step))
                    append (list noise (+ 0 (* i step)) g1))))))

(defun parse-data-line (str)
  (remove nil (mapcar (lambda (x) (read-from-string x nil nil))
                      (uiop:split-string str))))

(defun parse-data (file &key (start-line 0) (end-line nil))
  (with-open-file (in file)
    (loop :for line-number :from start-line
       :while (if end-line
                  (< line-number end-line)
                  t)
       :for line = (read-line in nil nil)
       :while line
       :collect (parse-data-line line))))

(defun write-list (filename list)
  (with-open-file (out filename :direction :output :if-exists :supersede :if-does-not-exist :create)
    (loop :for elem :in list
       :do (format out "~a~%" elem))))

(defun select-column (data column)
  (mapcar (lambda (x) (elt x column))
          data))

(defun select-columns (data start end)
  (loop :for i :from start :to end
     :collect (select-column data i)))

(defun transpose-lists (lists)
  (loop :for i :below (length (first lists))
     :collect (loop :for list :in lists
                 :collect (elt list i))))

(defun diff-data-line (data orig)
  (mapcar #'- data orig))

(defun diff-data-helper (list)
  (cons (caar list) (mapcar (lambda (x) (/ x (length list)))
                            (rest (reduce #'(lambda (x y) (mapcar #'- x y))
                                          list)))))

(defun avg-data (list)
  (cons (caar list) (mapcar (lambda (x) (/ x (length list)))
                            (rest (reduce #'(lambda (x y) (mapcar #'+ x y))
                                          list)))))
(defun avg-squared-data (list)
  (cons (caar list) (mapcar (lambda (x) (/ x (length list)))
                            (rest (reduce #'(lambda (x y) (mapcar (lambda (a b) (+ (expt a 2) (expt b 2))) x y))
                                          list)))))
(defun ms-data (list orig)
  (cons (caar list) (mapcar (lambda (x) (/ x (length list)))
                            (rest (reduce #'(lambda (x y) (mapcar #'+ x y))
                                          list)))))
(defun dir-equals-number (dir)
  (parse-integer (last-elt (pathname-directory dir)) :junk-allowed t))

(defun collect-numeric-directories (&optional (search-dir (uiop:getcwd)))
  (remove-if-not #'dir-equals-number
                 (uiop:subdirectories search-dir)))

(defun comment-line-p (line)
  "Checks whether line starts with #"
  (equal (elt line 0) #\#))

(defun get-time ()
  (multiple-value-bind (sec min hour) (get-decoded-time)
    (format nil "~2,'0d:~2,'0d:~2,'0d" hour min sec)))

(defun run (cmd)
  (format t "~&~a: ~a ~~$ ~a~%" (get-time) (uiop:getcwd) cmd)
  (let ((status (nth-value 2
                           (uiop:run-program cmd
                                             :output t
                                             :error-output t
                                             :ignore-error-status t))))
    (format t "~&Exit Code: ~a~%" status))
  (format t "~&~a: Done.~%" (get-time)))

(defvar *noise.gpl* "gnuplot -e \"set terminal png;
set pm3d map;
set output \\\"Noise.png\\\";
sp \\\"<gpo3 -re red_Noise.bin\\\" w pm3d;
print \\\"Noise done\\\";")

(defun generate-noise (parameter-file)
  (time (run (format nil "noise_gen ../~a ../inf_zero.bin" parameter-file)))
  (run "reduce_data Noise.bin")
                                        ;  (run *noise.gpl*)
  )

(defun link-to (path start end)
  (loop :for i :from start :to end
     :do (ensure-directories-exist (format nil "~a/" i))
       (if (uiop:absolute-pathname-p path)
           (run (format nil "ln -fs ~a/~a/Noise.bin ~:*~a/Noise.bin" path i))
           (run (format nil "ln -fs ../~a/~a/Noise.bin ~:*~a/Noise.bin" path i)))))

(defvar *rabi.gpl* "gnuplot -e \"set terminal png;
set output \\\"Rabi.png\\\";
p \\\"Rabi_1_0.txt\\\" u 1:2 w l, \\\"Rabi_1_0.txt\\\" u 1:3 w l;
print \\\"Rabi done\\\"\"")

(defun noise-solver (start end parameter-files)
  (loop :for param-file :in (ensure-list parameter-files)
     :do (loop :for dir :from start :to end
            :do
              (uiop:with-current-directory ((ensure-directories-exist
                                             (format nil "~a/~@[~a/~]"
                                                     dir
                                                     (pathname-name param-file))))
                (format t "Run in directory: ~a~%" (uiop:getcwd))
                (unless (every #'probe-file (list "../Noise.bin"))
                  (format t "Warning: Noise missing in ~a~%" dir))
                (format t "Solve~%")
                (time (run (format nil "noise_solver ../../~a" param-file)))))))

(defun chirp-solver (start end parameter-files)
  (loop :for param-file :in (ensure-list parameter-files)
     :do (loop :for dir :from start :to end
            :do
              (uiop:with-current-directory ((ensure-directories-exist
                                             (format nil "~a/~@[~a/~]"
                                                     dir
                                                     (pathname-name param-file))))
                (format t "Run in directory: ~a~%" (uiop:getcwd))
                (unless (every #'probe-file (list "../Noise.bin"))
                  (error "Noise missing in ~a~%" dir))
                (format t "Solve~%")
                (time (run "noise_solver more-chirps.xml"))))))

(defun diff-data (file1 file2)
  (let ((output-file (merge-pathnames (format nil "diff-~a.~a" (pathname-name file1) (pathname-type file1)) file1)))
    (with-open-file (in1 file1)
      (with-open-file (in2 file2)
        (with-open-file (out output-file
                             :direction :output :if-exists :supersede)
          (loop :for lines = (list (read-line in1 nil nil) (read-line in2 nil nil))
             :while (notany #'null lines)
             :do
               (if (comment-line-p (first lines))
                   (format out "~a~%" (first lines))
                   (format out "~{~a~t~}~%"
                           (diff-data-helper (mapcar #'parse-data-line lines))))))))
    output-file))

(defun average-data (start end filename subdirs)
  (loop :for subdir :in (ensure-list subdirs)
     :do (let ((filenames (loop :for num :from start :to end
                             :for path = (format nil "~a~a/~@[~a/~]~a" (uiop:getcwd) num subdir filename)
                             :collect path)))
           (let ((files (loop :for file :in filenames
                           :collect (open file :if-does-not-exist nil))))
             (unwind-protect
                  (progn (unless (every (lambda (fstream)
                                          (and (streamp fstream) (open-stream-p fstream)))
                                        files)
                           (error "Some files do not exist.~&~a" (loop :for file :in filenames

                                                                    :for s :in files
                                                                    :collect (cons file s))))
                         (with-open-file (out (ensure-directories-exist (format nil "avg-~a/~a" subdir filename))
                                              :direction :output :if-exists :supersede)
                           (loop :for lines = (loop :for file :in files
                                                 :collect (read-line file nil nil))
                              :while (notany #'null lines)
                              :do
                                (if (comment-line-p (first lines))
                                    (format out "~a~%" (first lines))
                                    (format out "~{~a~t~}~%"
                                            (avg-data (mapcar #'parse-data-line lines)))))
                           (format t "Average written into avg-~a/~a~%" subdir filename)))
               (mapcar (lambda (x) (when (and (streamp x) (open-stream-p x)) (close x)))
                       files))))))

(defun analyze-data (start end filename subdirs)
  (loop :for subdir :in (ensure-list subdirs)
     :do (let ((filenames (loop :for num :from start :to end
                             :for path = (format nil "~a~a/~@[~a/~]~a" (uiop:getcwd) num subdir filename)
                             :collect path)))
           (let ((files (loop :for file :in filenames
                           :collect (open file :if-does-not-exist nil))))
             (unwind-protect
                  (progn (unless (every (lambda (fstream)
                                          (and (streamp fstream) (open-stream-p fstream)))
                                        files)
                           (error "Some files do not exist.~&~a" (loop :for file :in filenames

                                                                    :for s :in files
                                                                    :collect (cons file s))))
                         (with-open-file (orig (format nil "fourier/~a/~a/" subdir filename))
                           (with-open-file (mean-out (ensure-directories-exist (format nil "ana-~a/mean-~a" subdir filename))
                                                     :direction :output :if-exists :supersede)
                             (with-open-file (var-out (ensure-directories-exist (format nil "ana-~a/var-~a" subdir filename))
                                                      :direction :output :if-exists :supersede)
                               (with-open-file (sigma-out (ensure-directories-exist (format nil "ana-~a/sigma-~a" subdir filename))
                                                          :direction :output :if-exists :supersede)
                                 (with-open-file (err-out (ensure-directories-exist (format nil "ana-~a/err-~a" subdir filename))
                                                          :direction :output :if-exists :supersede)
                                   (with-open-file (ms-out (ensure-directories-exist (format nil "ana-~a/ms-~a" subdir filename))
                                                           :direction :output :if-exists :supersede)
                                     (with-open-file (rms-out (ensure-directories-exist (format nil "ana-~a/rms-~a" subdir filename))
                                                              :direction :output :if-exists :supersede)
                                       (loop :for lines = (loop :for file :in files
                                                             :collect (read-line file nil nil))
                                          :for orig-line = (read-line orig nil nil)
                                          :while (and (notany #'null lines) (not (null orig-line)))
                                          :do
                                            (if (comment-line-p (first lines))
                                                (progn (format mean-out "~a~%" (first lines))
                                                       (format var-out "~a~%" (first lines))
                                                       (format sigma-out "~a~%" (first lines))
                                                       (format err-out "~a~%" (first lines))
                                                       (format ms-out "~a~%" (first lines))
                                                       (format rms-out "~a~%" (first lines)))
                                                (let* ((data-line (mapcar #'parse-data-line lines))
                                                       (orig-data (parse-data-line orig-line))
                                                       (avg-out (avg-data data-line))
                                                       (col0 (car avg-out))
                                                       (mean (cdr avg-out))
                                                       (mean-squared (cdr (avg-squared-data data-line)))
                                                       (variance (mapcar #'- mean-squared (mapcar (lambda (x) (expt x 2)) mean)))
                                                       (sigma (mapcar #'sqrt variance))
                                                       (diff-data (mapcar (lambda (x) (diff-data-line x orig-data)) data-line))
                                                       (err (cdr (avg-data diff-data)))
                                                       (ms (cdr (avg-squared-data diff-data)))
                                                       (rms (mapcar #'sqrt ms)))
                                                  (format mean-out "~{~a~t~}~%" (cons col0 mean))
                                                  (format var-out "~{~a~t~}~%" (cons col0 variance))
                                                  (format sigma-out "~{~a~t~}~%" (cons col0 sigma))
                                                  (format err-out "~{~a~t~}~%" (cons col0 err))
                                                  (format ms-out "~{~a~t~}~%" (cons col0 ms))
                                                  (format rms-out "~{~a~t~}~%" (cons col0 rms)))))))))))))
               (mapcar (lambda (x) (when (and (streamp x) (open-stream-p x)) (close x)))
                       files))))))

(defun average-binary-data (start end filename subdirs &key (average-function #'mean) (prefix "mean"))
  (loop :for subdir :in (ensure-list subdirs)
     :do (let ((filenames (loop :for num :from start :to end
                             :for path = (format nil "~a~a/~@[~a/~]~a" (uiop:getcwd) num subdir filename)
                             :collect path)))
           (let ((files (loop :for file :in filenames
                           :collect (open file :if-does-not-exist nil
                                          :element-type '(unsigned-byte 8)))))
             (unwind-protect
                  (progn (unless (every (lambda (fstream)
                                          (and (streamp fstream) (open-stream-p fstream)))
                                        files)
                           (error "Some files do not exist.~&~a" files))
                         (with-open-file (out (ensure-directories-exist (format nil "avg-~a/~@[~a-~]~a" subdir prefix filename))
                                              :direction :output :if-exists :supersede
                                              :element-type '(unsigned-byte 8))
                           (let ((header (read-header (first files))))
                             (loop :for file :in files
                                :do (file-position file (nself header)))
                             (write-bin out header)
                             (loop :for i :below (if (= (nDatatyp header) 0)
                                                     (* (nDimX header) (nDimY header) (nDimZ header))
                                                     (* (nDimX header) (nDimY header) (nDimZ header) 2))
                                :do (write-float64 out (funcall average-function (loop :for file :in files
                                                                                    :collect (read-float64 file))))))
                           (format t "Average written into avg-~a/~a~%" subdir filename)))
               (mapcar (lambda (x) (when (and (streamp x) (open-stream-p x)) (close x)))
                       files))))))

(defun abs-square (num)
  (expt (abs num) 2))

(defun collect-ensemble-files (start end filename &optional subdir)
  (loop :for i :from start :to end
     :collect (format nil "~d/~@[~a/~]~a" i subdir filename)))

(defun call-with-open-files (files fun &key (direction :input) (element-type 'base-char) if-exists if-does-not-exist more-args)
  "Opens every file in files and calls function fun with list of open streams as first argument and more-args as further arguments if any"
  (labels ((call-rec (files streams fun)
             (with-open-file (stream (first files)
                                     :direction direction :element-type element-type
                                     :if-exists if-exists :if-does-not-exist if-does-not-exist)
               (if (null (rest files))
                   (apply fun (cons (cons stream streams) more-args))
                   (call-rec (rest files) (cons stream streams) fun)))))
    (call-rec files nil fun)))

;; (defun average-binary-data2 (files output-file &key avg-fun)
;;   (flet ((avg-to-file-bin (streams output-file &key (avg-fun #'mean))
;;              (with-open-file (out (ensure-directories-exist output-file)
;;                                   :direction :output
;;                                   :if-exists :supersede
;;                                   :element-type '(unsigned-byte 8))
;;                (let ((header (read-header (first streams))))
;;                  (loop :for stream :in streams
;;                     :do (file-position stream (nself header)))
;;                  (write-header out header)
;;                  (loop :for i
;;                     :below (if (= (nDatatyp header) 0)
;;                                (* (nDimX header) (nDimY header) (nDimZ header))
;;                                (* (nDimX header) (nDimY header) (nDimZ header) 2))
;;                     :do (write-float64 out
;;                                        (apply avg-fun
;;                                               (loop :for file :in files
;;                                                  :collect (read-float64 file)))))
;;                  (format t "Average written into ~a~%" output-file)))))
;;     (call-with-open-files files avg-to-file-bin :element-type '(unsigned-byte 8))))

;; (defun average-binary-text2 (files output-file &key avg-fun)
;;   (flet ((avg-to-file-txt (streams output-file)
;;            (with-open-file (out (ensure-directories-exist output-file)
;;                                 :direction :output
;;                                 :if-exists :supersede)
;;              (loop :for lines = (loop :for file :in files
;;                                    :collect (read-line file nil nil))
;;                 :while (notany #'null lines)
;;                 :do
;;                   (if (comment-line-p (first lines))
;;                       (format out "~a~%" (first lines))
;;                       (format out "~{~a~t~}~%"
;;                               (avg-data (mapcar #'parse-data-line lines)))))
;;              (format t "Average written into ~a~%" output-file))))
;;     (call-with-open-files files avg-to-file-txt)))

(defmacro on-ensemble ((start end &optional (subdirs (list (list)))) &body body)
  (with-gensyms (subdir dir)
    `(loop :for ,subdir :in (ensure-list ,subdirs)
        :do (loop :for ,dir :from ,start :to ,end
               :do (uiop:with-current-directory ((ensure-directories-exist
                                                  (format nil "~a/~@[~a/~]"
                                                          ,dir
                                                          (when ,subdir
                                                            (pathname-name ,subdir)))))
                     ,@body)))))

(defun cmd (start end command subdirs)
  (assert (numberp start))
  (assert (numberp end))
  (on-ensemble (start end subdirs)
               (format t "~&  cmd:~a $ ~a~%" (uiop:getcwd) command)
               (run command)))

(defun generate-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (parameter-file (third args)))
        (loop :for dir :from start :to end
           :do (uiop:with-current-directory
                   ((ensure-directories-exist
                     (format nil
                             "~a/~@[~a/~]"
                             dir
                             nil)))
                 (format t "run in directory: ~a/~%" (uiop:getcwd))
                 (generate-noise (format nil "~a" parameter-file)))))
      (print "noise-solver generate <start> <end> <parameter-files*>")))

(defun filetype= (filename type)
  (assert (pathnamep filename))
  (assert (stringp type))
  (when (string= type (pathname-type (pathname filename)))
    filename))

(defun xml-file-p (filename)
  (filetype= filename "xml"))

(defun bin-file-p (filename)
  (filetype= filename "bin"))

(defun print-help ()
  (format t "~&noise-solver <param.xml> <start> <end> <mode> [<args>]

e.g.
noise-solver params.xml 1 20 generate
noise-solver params.xml 1 20 solve
noise-solver params.xml 1 20 average 100.000_1.bin gpo3 100.000_1.bin > 100.000_1.txt
noise-solver params.xml 1 20 cmd gpo3 100.000_1.bin > 100.000_1.txt
"))

(defun gather (start end filename &optional subdir )
  (ensure-directories-exist "gather/")
  (loop :for dir :from start :to end
     :with file = (format nil "~a/~@[~a/~]~a" dir subdir filename)
     :do (uiop:copy-file file
                         (format nil "~a/~a" "gather" (tag-pathname filename dir)))))

(defun tag-pathname (pathname &rest tags)
  (merge-pathnames (format nil "~a~{-~a~}" (pathname-name pathname) tags)
                   (pathname pathname)))

(defun read-number (in)
  (let* ((pos (file-position in))
         (line (read-line in nil nil)))
    (when line
      (multiple-value-bind (num end) (parse-float:parse-float line :junk-allowed t)
        (when (/= (+ pos end 1) (file-position in))
          (file-position in (+ pos end)))
        num))))

(defun rms (sample)
  "Calculates root means squared of sample"
  (if (null sample)
      0
      (sqrt (/ (reduce (lambda (base elem)
                         (+ base (expt elem 2)))
                       sample :initial-value 0)
               (length sample)))))

(defun diff-visibility (file1 file2)
  (with-open-file (out (format nil "diff-~a.txt" (pathname-name file2))
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
    (with-open-file (in1 file1)
      (with-open-file (in2 file2)
        (loop
           :for line2 = (read-line in2 nil nil)
           :while line2
           :for (x2 y2) = (parse-data-line line2)
           :do (loop
                  :for line1 = (read-line in1 nil nil)
                  :while line1
                  :for (x1 y1) = (parse-data-line line1)
                  :when (= x1 x2)
                  :unless (> 65000 x1 45000)
                  :do
                    (format out "~a ~a~%"
                            x1
                            (rms y1))
                    (return-from nil))
             (file-position in1 0))))))

(defun diff-visibility2 (file1 file2)
  (with-open-file (out (format nil "diff2-~a.txt" (pathname-name file2))
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
    (with-open-file (in1 file1)
      (with-open-file (in2 file2)
        (loop
           :for line2 = (read-line in2 nil nil)
           :while line2
           :for (x2 . y2) = (parse-data-line line2)
           :do (loop
                  :for line1 = (read-line in1 nil nil)
                  :while line1
                  :for (x1 y1) = (parse-data-line line1)
                  :when (= x1 x2)
                  :unless (> 65000 x1 45000)
                  :do
                    (let ((diff (mapcar (lambda (y) (abs (- y y1)))
                                        y2)))
                      (format out "~a ~a ~a ~a ~a~%"
                              x1
                              (rms diff)
                              (mean diff)
                              (variance diff)
                              (standard-deviation diff))
                      (return-from nil)))
             (file-position in1 0))))))

(defun visibility (path)
  (uiop:with-current-directory (path)
    (run (format nil "interpol_chirp Chirp_5.txt"))
    (with-open-file (in "Chirp_interpol.txt")
      (let* ((data (loop :for num = (read-number in)
                      :while num
                      :collect (read-number in)))
             (imax (apply #'max data))
             (imin (apply #'min data)))
        (/ (- imax imin) (+ imax imin))))))

(defun get-bragg-visibility (strength duration g1)
  (with-open-file (out (format nil "visibility-fourier-g1-~a.txt" g1) :direction :output
                       :if-exists :supersede :if-does-not-exist :create)
    (loop :for time :in (ensure-list duration)
       :for file = (format nil "fourier/if-~f-chirp-~d~@[-g1-~f~]/Chirp_5.txt"
                           strength time g1)
       :do (format out "~d ~f~%" time (visibility file)))))

(defun get-visibility (strength &optional (duration (list 1000 2000 3000 5000 10000 14000 21000)) (g1 nil))
  (with-open-file (out (format nil "visibility-~a~@[-~a~].txt" strength g1)
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
    (loop :for time :in (ensure-list duration)
       :do (format out "~a ~a~%" time (visibility (format nil "avg-if-~@[~a-~]chirp-~a-~@[g1-~a-~]Chirp_5.txt" strength time g1))))))

(defun get-all-visibility (start end strength &optional (duration (list 1000 2000 3000 5000 10000 14000 21000)) (g1 nil))
  (with-open-file (out (format nil "visibility-all-~a~@[-~a~].txt" strength g1)
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
    (loop :for time :in (ensure-list duration)
       :do (format out "~a~{ ~a~}~%" time (loop :for num :from start :to end
                                             :collect (visibility (format nil "~a/if-~@[~f-~]chirp-~d~@[-g1-~f~]/Chirp_5.txt" num strength time g1)))))))

(defun get-avg-visibility (start end strength &optional (duration (list 1000 2000 3000 5000 10000 14000 21000)) (g1 nil))
  (with-open-file (out (format nil "visibility-avg-~a~@[-~a~].txt" strength g1)
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
    (loop :for time :in (ensure-list duration)
       :do (format out "~a ~a~%" time (mean (loop :for num :from start :to end
                                               :collect (visibility (format nil "~a/if-~@[~f-~]chirp-~d~@[-g1-~f~]/Chirp_5.txt" num strength time g1))))))))


(defun average-files (start end filename strength
                      &optional
                        (duration (list 1000 2000 3000 5000 10000 14000 21000)) (g1 nil))
  (loop :for time :in duration
     :for subdir = (format nil "if-~@[~a-~]chirp-~a~@[-g1-~a~]" strength time g1)
     :do (average-data start end filename subdir)))

(defvar *duration* (list 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000 65000 70000 75000 80000))
(defvar *duration2* (list 10000 15000 20000 25000 30000 35000 40000 45000 50000 80000))

(defun generate-xml (strength duration &optional (chirp t) (g1 0.0))
  (loop :for str :in (ensure-list strength)
     :do (loop :for time :in (ensure-list duration)
            :do (loop :for g :in (ensure-list g1)
                   :do (with-open-file (out (format nil "if-~@[~f-~]~@[chirp-~]~*~d~@[-g1-~f~].xml"
                                                    str chirp time g)
                                            :direction :output
                                            :if-exists :supersede
                                            :if-does-not-exist :create)
                         (format out "<SIMULATION>
  <DIM>1</DIM>
  <FILENAME>../../inf_0.000_1.bin</FILENAME>
  <FILENAME_2>../../inf_zero.bin</FILENAME_2>
  <FILENAME_3>../../inf_0.000_1.bin</FILENAME_3>
  <FILENAME_4>../../inf_zero.bin</FILENAME_4>
  <NOISE>../Noise.bin</NOISE>
  <CONSTANTS>
    <laser_k>8.05289</laser_k>
    <laser_k_2>8.05289</laser_k_2>
    <laser_domh>0.0471239</laser_domh>
    <laser_domh_2>0.0471239</laser_domh_2>
    <laser_dk>0</laser_dk>
    <rabi_threshold>4</rabi_threshold>
    <Noise_Amplitude>~f</Noise_Amplitude>
  </CONSTANTS>
  <VCONSTANTS>
    <Amp_1>-14.0496,-14.0496</Amp_1>
    <Amp_2>-14.0496,-14.0496</Amp_2>
    <Alpha_1>0.000365368,0.000365368,0.000365368</Alpha_1>
    <Alpha_2>0.000365368,0.000365368,0.000365368</Alpha_2>
    <Delta_L>0,-6283.19,0,0</Delta_L>
    <GS_1>~f,~:*~f,0,0,0,0</GS_1>
    <GS_2>~:*~f,~:*~f,0,0,0,0</GS_2>
    <GS_3>0,0,0,0,0,0</GS_3>
    <GS_4>0,0,0,0,0,0</GS_4>
    <GS_5>0,0,0,0,0,0</GS_5>
    <GS_6>0,0,0,0,0,0</GS_6>
    <Beta>0.0,0.0,0.0</Beta>
  </VCONSTANTS>
  <ALGORITHM>
    <NX>16384</NX>
    <NY>1</NY>
    <NZ>1</NZ>
    <XMIN>-320</XMIN>
    <XMAX>320</XMAX>
    <NK>25</NK>
    <NA>700</NA>
    <EPSILON>1e-6</EPSILON>
  </ALGORITHM>
  <SEQUENCE>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"packed\" pn_freq=\"last\" rabi_output_freq=\"each\">100</bragg_ad>
    <freeprop dt=\"1\" Nk=\"100\" output_freq=\"packed\" pn_freq=\"last\">~d</freeprop>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"packed\" pn_freq=\"last\" rabi_output_freq=\"last\" >200</bragg_ad>
    <freeprop dt=\"1\" Nk=\"100\" output_freq=\"packed\" pn_freq=\"last\">~:*~d</freeprop>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"packed\" pn_freq=\"last\" rabi_output_freq=\"each\"~@[ chirp_mode=\"1\" no_of_chirps=\"10\"~]>100</bragg_ad>
  </SEQUENCE>
</SIMULATION>"
                                 str g (/ time 2) chirp))))))

(defun generate-fast-xml (strength duration &optional (chirp t) (g1 0.0) (dir 45000))
  (loop :for str :in (ensure-list strength)
     :do (loop :for time :in (ensure-list duration)
            :do (loop :for g :in (ensure-list g1)
                   :do (with-open-file (out (format nil "if-~@[~f-~]~@[chirp-~]~*~d~@[-g1-~f~].xml"
                                                    str chirp time g)
                                            :direction :output
                                            :if-exists :supersede
                                            :if-does-not-exist :create)
                         (format out "<SIMULATION>
  <DIM>1</DIM>
  <FILENAME>../if-~f-chirp-~d-g1-~f/~d.000_1.bin</FILENAME>
  <FILENAME_2>../../inf_zero.bin</FILENAME_2>
  <FILENAME_3>../~4:*if-~f-chirp-~d-g1-~f/~d.000_1.bin</FILENAME_3>
  <FILENAME_4>../../inf_zero.bin</FILENAME_4>
  <NOISE>../Noise.bin</NOISE>
  <CONSTANTS>
    <laser_k>8.05289</laser_k>
    <laser_k_2>8.05289</laser_k_2>
    <laser_domh>0.0471239</laser_domh>
    <laser_domh_2>0.0471239</laser_domh_2>
    <laser_dk>0</laser_dk>
    <rabi_threshold>4</rabi_threshold>
    <Noise_Amplitude>~f</Noise_Amplitude>
  </CONSTANTS>
  <VCONSTANTS>
    <Amp_1>-14.0496,-14.0496</Amp_1>
    <Amp_2>-14.0496,-14.0496</Amp_2>
    <Alpha_1>0.000365368,0.000365368,0.000365368</Alpha_1>
    <Alpha_2>0.000365368,0.000365368,0.000365368</Alpha_2>
    <Delta_L>0,-6283.19,0,0</Delta_L>
    <GS_1>~f,~:*~f,0,0,0,0</GS_1>
    <GS_2>~:*~f,~:*~f,0,0,0,0</GS_2>
    <GS_3>0,0,0,0,0,0</GS_3>
    <GS_4>0,0,0,0,0,0</GS_4>
    <GS_5>0,0,0,0,0,0</GS_5>
    <GS_6>0,0,0,0,0,0</GS_6>
    <Beta>0.0,0.0,0.0</Beta>
  </VCONSTANTS>
  <ALGORITHM>
    <NX>16384</NX>
    <NY>1</NY>
    <NZ>1</NZ>
    <XMIN>-320</XMIN>
    <XMAX>320</XMAX>
    <NK>25</NK>
    <NA>700</NA>
    <EPSILON>1e-6</EPSILON>
  </ALGORITHM>
  <SEQUENCE>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"packed\" pn_freq=\"last\" rabi_output_freq=\"last\" >200</bragg_ad>
    <freeprop dt=\"1\" Nk=\"100\" output_freq=\"packed\" pn_freq=\"last\">~d</freeprop>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"packed\" pn_freq=\"last\" rabi_output_freq=\"each\"~@[ chirp_mode=\"1\" no_of_chirps=\"10\"~]>100</bragg_ad>
  </SEQUENCE>
</SIMULATION>"
                                 str dir g (+ 100 (/ time 2)) str g (/ time 2) chirp))))))

(defun generate-more-chirp-xml (start end strength duration &optional (chirp 30) (g1 0.0))
  (loop :for str :in (ensure-list strength)
     :do (loop :for time :in (ensure-list duration)
            :do (loop :for g :in (ensure-list g1)
                   :do (loop :for run :from start :to end
                          :do (with-open-file (out (format nil "~a/if-~@[~f-~]~@[chirp-~]~*~d~@[-g1-~f~]/more-chirps.xml"
                                                           run str chirp time g)
                                                   :direction :output
                                                   :if-exists :supersede
                                                   :if-does-not-exist :create)
                                (format out "<SIMULATION>
  <DIM>1</DIM>
  <FILENAME>~a.000_1.bin</FILENAME>
  <FILENAME_2>../../inf_zero.bin</FILENAME_2>
  <FILENAME_3>~:*~a.000_1.bin</FILENAME_3>
  <FILENAME_4>../../inf_zero.bin</FILENAME_4>
  <NOISE>../Noise.bin</NOISE>
  <CONSTANTS>
    <laser_k>8.05289</laser_k>
    <laser_k_2>8.05289</laser_k_2>
    <laser_domh>0.0471239</laser_domh>
    <laser_domh_2>0.0471239</laser_domh_2>
    <laser_dk>0</laser_dk>
    <rabi_threshold>4</rabi_threshold>
    <Noise_Amplitude>~f</Noise_Amplitude>
  </CONSTANTS>
  <VCONSTANTS>
    <Amp_1>-14.0496,-14.0496</Amp_1>
    <Amp_2>-14.0496,-14.0496</Amp_2>
    <Alpha_1>0.000365368,0.000365368,0.000365368</Alpha_1>
    <Alpha_2>0.000365368,0.000365368,0.000365368</Alpha_2>
    <Delta_L>0,-6283.19,0,0</Delta_L>
    <GS_1>~f,~:*~f,0,0,0,0</GS_1>
    <GS_2>~:*~f,~:*~f,0,0,0,0</GS_2>
    <GS_3>0,0,0,0,0,0</GS_3>
    <GS_4>0,0,0,0,0,0</GS_4>
    <GS_5>0,0,0,0,0,0</GS_5>
    <GS_6>0,0,0,0,0,0</GS_6>
    <Beta>0.0,0.0,0.0</Beta>
  </VCONSTANTS>
  <ALGORITHM>
    <NX>16384</NX>
    <NY>1</NY>
    <NZ>1</NZ>
    <XMIN>-320</XMIN>
    <XMAX>320</XMAX>
    <NK>25</NK>
    <NA>700</NA>
    <EPSILON>1e-6</EPSILON>
  </ALGORITHM>
  <SEQUENCE>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"each\" pn_freq=\"last\" rabi_output_freq=\"each\" chirp_mode=\"1\" no_of_chirps=\"~a\">100</bragg_ad>
  </SEQUENCE>
</SIMULATION>"
                                        (+ time 300) str g chirp)))))))

(defun bragg-solve (param-files)
  (loop :for file :in (ensure-list param-files)
     :do (uiop:with-current-directory
             ((ensure-directories-exist
               (format nil "fourier/~a/" (pathname-name file))))
           (run (concatenate 'string "bragg ../../" file)))))

(defun bragg-solve-helper (args)
  (if args
      (bragg-solve args)
      (print "noise-solver bragg-solve <param-files*>")))

(defun solve-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (param-files (subseq args 2)))
        (noise-solver start end param-files))
      (print "noise-solver solve <start> <end> <param-files*>")))

(defun chirp-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (param-files (subseq args 2)))
        (chirp-solver start end param-files))
      (print "noise-solver solve <start> <end> <param-files*>")))

(defun average-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (filename (third args))
            (subdirs (mapcar #'pathname-name (subseq args 3))))
        (if (string= (pathname-type filename) "bin")
            (average-binary-data start end filename subdirs )
            (average-data start end filename subdirs)))
      (print "noise-solver average <start> <end> <filename> <subdirs*>")))

(defun cmd-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (command (third args))
            (subdirs (mapcar #'pathname-name (subseq args 3))))
        (cmd start end command subdirs))
      (print "noise-solver cmd <start> <end> <command> <subdirs*>")))

(defun fetch (start end filename subdirs)
  (ensure-directories-exist "fetched/")
  (loop :for num :from start :to end
     :do (loop :for subdir :in (ensure-list subdirs)
            :do (uiop:copy-file (format nil "~a/~a/~a" num subdir filename)
                                (format nil "fetched/~a-~a-~a" num subdir filename)))))

(defun fetch-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (filename (third args))
            (subdirs (mapcar #'pathname-name (subseq args 3))))
        (fetch start end filename subdirs))
      (print "noise-solver fetch <start> <end> <filename> <subdirs*>")))


(defun main (&optional (argv (uiop:command-line-arguments)))
  (when (null argv)
    (format t "~&No commands given.~%")
    (uiop:quit))
  (switch ((first argv) :test #'string=)
          ("solve" (solve-helper (rest argv)))
          ("bragg-solve" (bragg-solve-helper (rest argv)))
          ("chirp-solve" (chirp-helper (rest argv)))
          ("generate" (generate-helper (rest argv)))
          ("average" (average-helper (rest argv)))
          ("cmd" (cmd-helper (rest argv)))
          ("fetch" (fetch-helper (rest argv)))
          ("repl" (format t "~%> ") (print (eval (read))))
          (t (format t "Command not known: ~a~%" (first argv))
             (uiop:quit))))

;; (defun ana-data ((data atus-vec))
;;   (assert (= (nDims data) 1))
;;   (loop
;;      :for el :across (data data)
;;      :for i upfrom 0
;;      :for x = (+ (xMin data) (* (dx data) i))
;;      :with mean = 0.0
;;      :with variance = 0.0
;;      :with sum = 0.0
;;      :do
;;        (incf mean (* x el))
;;        (incf variance (* x x el))
;;        (incf sum el)
;;      :finally  (return (list (/ mean sum) (/ variance sum) sum))))

(defun set-pathtype (file type)
  (make-pathname :name (pathname-name file) :directory (pathname-directory file) :type type))

(defun gpo3-on-dir (dir)
  (let ((files (remove-if-not (lambda (file) (string= (pathname-type file) "bin"))
		              (uiop:directory-files dir))))
    (mapcar (lambda (file) (uiop:run-program (format nil "/home/strelox/bin/gpo3 ~a > ~a" file (set-pathtype file "txt"))))
            files)
    t))
