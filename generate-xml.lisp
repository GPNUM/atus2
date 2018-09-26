#!/usr/bin/sbcl --script

(setf *read-default-float-format* 'double-float)

(defun safe-read ()
  "Somewhat hardened against malicious read-macros"
  (with-standard-io-syntax
    (let* ((*read-eval* nil))
      (read))))

(defun ask-for-number (prompt)
  (let ((entry))
    (loop
       (fresh-line)
       (format t prompt)
       (finish-output)
       (setf entry (safe-read))
       (when (numberp entry)
         (return entry))
       (format t "~&Not a valid number: ~a" entry))))

(defparameter seq-list (list "free propagation"
                             "interferometer"))

(let ((xml-form "<SIMULATION>
  <DIM>~d</DIM>
  <FILENAME>0.0000_1.bin</FILENAME>
  <FILENAME_2>0.0000_2.bin</FILENAME_2>
  <CONSTANTS>
    <Beta>-0.00134248</Beta>
    <laser_k>8.05289</laser_k>
    <laser_k_2>8.05289</laser_k_2>
    <laser_domh>0.0471239</laser_domh>
    <laser_domh_2>0.0471239</laser_domh_2>
    <laser_dk>0</laser_dk>
    <rabi_threshold>4</rabi_threshold>
    <chirp>0</chirp>
  </CONSTANTS>
  <VCONSTANTS>
    <Amp_1>-14.0496,-14.0496</Amp_1>
    <Amp_2>-14.0496,-14.0496</Amp_2>
    <Alpha_1>~,6e~{,~,6e~}</Alpha_1>~:[~1*~;~:*~%    <Alpha_2>~,6e~{,~,6e~}</Alpha_2>~]
    <Delta_L>0,-6283.19</Delta_L>
    <GS_1>~,6e,~,6e</GS_1>
    <GS_2>~,6e,~,6e</GS_2>
    <Beta>~,6e~{,~,6e~}</Beta>~:[~1*~;~:*~%    <Beta_2>~,6e~{,~,6e~}</Beta_2>~]
  </VCONSTANTS>
  <ALGORITHM>
    <NK>25</NK>
    <NA>700</NA>
  </ALGORITHM>
  ~-[<SEQUENCE>
    <freeprop dt=\"0.2\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\">1000</freeprop>
  </SEQUENCE>~;<SEQUENCE>
    <bragg_ad dt=\"0.1\" Nk=\"100\" output_freq=\"last\" pn_freq=\"last\" rabi_output_freq=\"each\" >100</bragg_ad>
    <freeprop dt=\"0.1\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\">7000</freeprop>
    <bragg_ad dt=\"0.1\" Nk=\"100\" output_freq=\"last\" pn_freq=\"last\" rabi_output_freq=\"last\" >200</bragg_ad>
    <freeprop dt=\"0.1\" Nk=\"10\" output_freq=\"last\" pn_freq=\"last\">7000</freeprop>
    <bragg_ad dt=\"0.1\" Nk=\"100\" output_freq=\"last\" pn_freq=\"last\" rabi_output_freq=\"each\" phase_scan=\"20\">100</bragg_ad>
  </SEQUENCE>~]
</SIMULATION>
"))


  (defun write-xml (filename dimension alpha alpha2 beta beta2 gs sequence)
    (with-open-file (out filename :direction :output :if-exists :supersede :if-does-not-exist :create)
      (format out xml-form dimension (car alpha) (cdr alpha) (car alpha2) (cdr alpha2) (nth 0 gs) (nth 1 gs) (nth 2 gs) (nth 3 gs) (car beta) (cdr beta) (car beta2) (cdr beta2) sequence)
      (format t "~&~a written.~%" filename))))


(defun generate-xml ()
  (let ((hbar 1.0545718e-34) ; reduced Planck constant
        (kg/u 1.66053904e-27) ; u -> kg conversion constant
        (two-species-p nil)
        (dimension)
        (length 1e-6) ; Length scale
        (time 1e-6) ; Time scale
        (gravity 0)
        (mass1 87) ; Mass of species 1
        (mass2 85) ; Mass of species 2
        (alpha) ; Kinetic factor in Hamiltonian
        (alpha2) ; Kinetic factor in Hamiltonian for second species
        (beta) ; Gravity
        (beta2)
        (gs) ; Non-linearity
        (sequence))

    ;; local functions
    (flet ((calculate-alpha (mass)
             (/ (* hbar time ) 2 (* mass kg/u) (expt length 2)))
           (calculate-beta (mass)
             (/ (* time length gravity (* mass kg/u))
                hbar))
           (calculate-gs (mass1 mass2 a)
             (/ (* 2 pi hbar time a (+ mass1 mass2))
                mass1
                mass2
                (expt length dimension)))
           (print-settings (dim)
             (format t "~&~%Settings along dimension: ~a" (nth dim (list "x" "y" "z")))
             (format t "
Current Setting:
  Length scale (in m): ~a
  Time scale (in s): ~a
  Gravity (in m/s^2): ~a
  ~:[Mass (in u): ~a~1*~;Mass1 (in u): ~a~%  Mass2 (in u): ~a~]
What would you change?
1 - Length
2 - Time
3 - Gravity
~:[4 - Mass~;4 - Mass1~%5 - Mass2~]
p - print current settings
c - continue
q - abort~%" length time gravity two-species-p mass1 mass2 two-species-p)))

      ;; Start
      (setf two-species-p (yes-or-no-p "Two species?"))

      ;; Set number of dimensions
      (loop
         (format t "Number of dimensions? (1-3) ")
         (finish-output)
         (setf dimension (parse-integer (read-line) :junk-allowed t))
         (when (and (numberp dimension)
                    (some (lambda (x) (= dimension x)) (list 1 2 3)))
           (return))
         (format t "~&Dimension not supported: ~a~%" dimension))

      ;; Set configuration per dimension
      (loop :for dim :from 0 :below dimension
         :do
           (print-settings dim)
           (format t "> ")
           (finish-output)
           (loop :for c = (read-line)
              :until (string= c "c")
              :do
                (cond ((string= c "1")
                       (setf length (ask-for-number "Length scale (in m): "))
                       (format t "~&Length scale updated.~%"))
                      ((string= c "2")
                       (setf time (ask-for-number "Time scale (in s): "))
                       (format t "~&Time scale updated.~%"))
                      ((string= c "3")
                       (setf gravity (ask-for-number "Gravity (in m/s^2): "))
                       (format t "~&Gravity updated.~%"))
                      ((string= c "4")
                       (setf mass1 (ask-for-number "Mass (in atomic mass units, u): "))
                       (unless two-species-p
                         (setf mass2 mass1))
                       (format t "~&Mass updated.~%"))
                      ((and two-species-p (string= c "5"))
                       (setf mass2 (ask-for-number "Mass2 (in atomic mass units, u): "))
                       (format t "~&Mass2 updated.~%"))
                      ((string= c "p")
                       (print-settings dim))
                      ((string= c "q")
                       (format t "~&Quit.~%")
                       (quit))
                      (t (format t "~&Unknown command.~%")))
                (format t "> ")
                (finish-output))
           (setf alpha (cons (calculate-alpha mass1) alpha))
           (setf beta (cons (calculate-beta mass1) beta))
           (when two-species-p
             (setf alpha2 (cons (calculate-alpha mass2) alpha2))
             (setf beta2 (cons (calculate-beta mass2) beta2))))

      ;; Set nonlinearity
      (setf gs (if (yes-or-no-p "Calculate nonlinearity?")
                   (list (calculate-gs mass1 mass1 (ask-for-number "a_11: "))
                         (calculate-gs mass1 mass2 (ask-for-number "a_12: "))
                         (calculate-gs mass2 mass1 (ask-for-number "a_21: "))
                         (calculate-gs mass2 mass2 (ask-for-number "a_22: ")))
                   (list 0.0 0.0 0.0 0.0)))

      ;; Set sequence
      (format t "~&Select Sequence:~%~{~a - ~a~%~}> "
              (loop :for seq :in seq-list
                 :for i :from 0
                 :append (list i seq)))
      (finish-output)
      (setf sequence (loop :for c = (parse-integer (read-line) :junk-allowed t)
                        :do
                          (if (and (numberp c) (>= c 0) (< c (length seq-list)))
                              (return c)
                              (progn (format t "~&Not valid: ~a Try again.~%> " c)
                                     (finish-output)))))

      (write-xml "output.xml" dimension alpha alpha2 beta beta2 gs sequence))))

;; Call program
(generate-xml)
