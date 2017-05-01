;;;; noise-solver.lisp

(in-package #:noise-solver)

(setf *read-default-float-format* 'double-float)

(defun parse-data-line (str)
  (remove nil (mapcar (lambda (x) (read-from-string x nil nil))
                      (uiop:split-string str))))

(defun avg-data (list)
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
  (multiple-value-bind (sec min hour) (decode-universal-time (get-universal-time))
    (format nil "~a:~a:~a" hour min sec)))

(defun run (cmd)
  (format t "~&~a ~~$ ~a~%" (uiop:getcwd) cmd)
  (uiop:run-program cmd
                    :output nil
                    :error-output nil))

(defvar *noise.gpl* "gnuplot -e \"set terminal png;
set pm3d map;
set output \\\"Noise.png\\\";
sp \\\"<gpo3 -re red_Noise.bin\\\" w pm3d;
print \\\"Noise done\\\";")

(defun generate-noise (parameter-file)
  (run (format nil "noise_gen ../~a ../inf_zero.bin" parameter-file))
  (run "reduce_data Noise.bin")
;  (run *noise.gpl*)
  )

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
                (error "Noise missing in ~a~%" dir))
              (format t "Solve~%")
              (time (run (format nil "noise_solver ../../~a" param-file)))))))

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
                         (with-open-file (out (format nil "avg-~@[~a-~]~a" subdir filename)
                                              :direction :output :if-exists :supersede)
                           (loop :for lines = (loop :for file :in files
                                                 :collect (read-line file nil nil))
                              :while (notany #'null lines)
                              :do
                              (if (comment-line-p (first lines))
                                  (format out "~a~%" (first lines))
                                  (format out "~{~a~t~}~%"
                                          (avg-data (mapcar #'parse-data-line lines)))))
                           (format t "Average written into avg-~a-~a~%" subdir filename)))
               (mapcar (lambda (x) (when (and (streamp x) (open-stream-p x)) (close x)))
                       files))))))

(defun average-binary-data (start end filename subdirs)
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
                         (with-open-file (out (format nil "avg-~@[~a-~]~a" subdir filename)
                                              :direction :output :if-exists :supersede
                                              :element-type '(unsigned-byte 8))
                           (let ((header (read-header (first files))))
                             (loop :for file :in files
                                :do (file-position file (nself header)))
                             (write-header out header)

                             (loop :for i :below (if (= (nDatatyp header) 0)
                                                     (* (nDimX header) (nDimY header) (nDimZ header))
                                                     (* (nDimX header) (nDimY header) (nDimZ header) 2))
                                 :do (write-float64 out (mean (loop :for file :in files
                                                           :collect (read-float64 file))))))
                           (format t "Average written into avg-~a-~a~%" subdir filename)))
               (mapcar (lambda (x) (when (and (streamp x) (open-stream-p x)) (close x)))
                       files))))))

(defun cmd (start end command subdirs)
  (assert (numberp start))
  (assert (numberp end))
  (loop :for subdir :in (ensure-list subdirs)
     :do (loop :for dir :from start :to end
            :do (uiop:with-current-directory ((format nil "~a/~a/" dir subdir))
                  (format t "~&  cmd:~a $ ~a~%" (uiop:getcwd) command)
                  (run command)))))

(defun generate-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (parameter-files (subseq args 2)))
        (loop :for parameter-file :in parameter-files
              :do (loop :for dir :from start :to end
                        :do
                           (uiop:with-current-directory
                               ((ensure-directories-exist
                                 (format nil "~a/" dir)))
                             (format t "Run in directory: ~a/~%" (uiop:getcwd))
                             (generate-noise (format nil "~a" parameter-file))))))
      (print "noise-solver generate <start> <end> <parameter-files*>")))

(defun xml-file-p (filename)
  (when (string= "xml" (pathname-type (pathname filename)))
    filename))

(defun bin-file-p (filename)
  (when (string= "bin" (pathname-type (pathname filename)))
    filename))

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
  (let ((num (read in nil nil)))
    (if (or (numberp num) (null num))
        num
        (error "Not a number: ~a" num))))

(defun read-number-from-string (string &key (eof-error-p t) (start 0))
  (let ((num (read-from-string string eof-error-p nil :start start)))
    (if (numberp num)
        num
        (error "Not a number: ~a" num))))

(defun visibility (path)
  (with-open-file (in path)
    (read-line in) ;; Skip comment
    (let* ((data (loop :for num = (read-number in)
                    :while num
                    :do (read-number in)
                    :collect (read-number in)
                    :do (read-number in)
                    (read-number in)))
           (imax (apply #'max data))
           (imin (apply #'min data)))
      (/ (- imax imin) (+ imax imin)))))

(defun get-bragg-visibility (strength duration g1)
  (with-open-file (out "visibility-fourier-g1-0.0001.txt" :direction :output
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

(defun average-files (start end filename strength
                      &optional
                        (duration (list 1000 2000 3000 5000 10000 14000 21000)) (g1 nil))
  (loop :for time :in duration
     :for subdir = (format nil "if-~@[~a-~]chirp-~a~@[-g1-~a~]" strength time g1)
     :do (average-data start end filename subdir)))

(defun generate-xml (strength duration &optional (chirp t) (g1 0.0) (noise 2.0))
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
  <DX_NOISE>../dxNoise.bin</DX_NOISE>
  <DX2_NOISE>../dx2Noise.bin</DX2_NOISE>
  <CONSTANTS>
    <laser_k>8.05289</laser_k>
    <laser_k_2>8.05289</laser_k_2>
    <laser_domh>0.0471239</laser_domh>
    <laser_domh_2>0.0471239</laser_domh_2>
    <laser_dk>0</laser_dk>
    <rabi_threshold>4</rabi_threshold>
    <Noise_Amplitude>~f</Noise_Amplitude>
    <Noise_xCorr>~f</Noise_xCorr>
    <Noise_tCorr>~:*~f</Noise_tCorr>
    <Noise_Sigma>1.0</Noise_Sigma>
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
    <NX>8192</NX>
    <NY>1</NY>
    <NZ>1</NZ>
    <XMIN>-160</XMIN>
    <XMAX>160</XMAX>
    <NK>25</NK>
    <NA>700</NA>
    <EPSILON>1e-6</EPSILON>
  </ALGORITHM>
  <SEQUENCE>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\" rabi_output_freq=\"each\">100</bragg_ad>
    <freeprop dt=\"0.1\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\">~d</freeprop>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\" rabi_output_freq=\"last\" >200</bragg_ad>
    <freeprop dt=\"0.1\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\">~:*~d</freeprop>
    <bragg_ad dt=\"0.1\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\" rabi_output_freq=\"each\"~@[ chirp_mode=\"1\" no_of_chirps=\"10\"~]>100</bragg_ad>
  </SEQUENCE>
</SIMULATION>"
                                 str noise g (/ time 2) chirp))))))

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

(defun average-helper (args)
  (if args
      (let ((start (parse-integer (first args)))
            (end (parse-integer (second args)))
            (filename (third args))
            (subdirs (mapcar #'pathname-name (subseq args 3))))
        (if (string= (pathname-type filename) "bin")
            (average-binary-data start end filename subdirs)
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
    (error "No commands given."))
  (switch ((first argv) :test #'string=)
    ("solve" (solve-helper (rest argv)))
    ("bragg-solve" (bragg-solve-helper (rest argv)))
    ("generate" (generate-helper (rest argv)))
    ("average" (average-helper (rest argv)))
    ("cmd" (cmd-helper (rest argv)))
    ("fetch" (fetch-helper (rest argv)))
    (t (error "Command not known: ~a" (first argv)))))
