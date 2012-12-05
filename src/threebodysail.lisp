(in-package :bld-threebodysail)

(defparameter *double-eps* 1d-10)

(defun is= (a b)
  (< (abs (- (abs a) (abs b))) *double-eps*))

(defun is> (a b)
  (> (- a b) *double-eps*))

(defun is>= (a b)
  (or (is> a b) (is= a b)))

(defmethod coincident ((a ve3) (b ve3))
  (is= (norme2 (*x2 a b)) 0d0))

(defmethod codirectional ((a ve3) (b ve3))
  (and (coincident a b)
       (is> (scalar (*i a b)) 0d0)))

(defmethod dot ((a ve3) (b ve3))
  (scalar (*i a b)))

(defun get-window-axis-properties (window axis)
  (cond ((eq axis :x) (slot-value (slot-value window 'cl-plplot::x-axis) 'cl-plplot::properties))
	((eq axis :y) (slot-value (slot-value window 'cl-plplot::y-axis) 'cl-plplot::properties))
	(t (warn "Axis must be one of :x or :y") nil)))

(defun gradv (res r12 mu1 mu2)
  (let* ((mu (/ mu2 (+ mu1 mu2)))
	 (re1 (ve3 :c1 (* -1 r12 mu)))
	 (re2 (ve3 :c1 (* r12 (- 1 mu))))
	 (r1s (- res re1))
	 (r2s (- res re2))
	 (r1ss (norme r1s))
	 (r2ss (norme r2s))
	 (a1g (if (is> r1ss 0)
		  (* (- (/ mu1 (expt r1ss 3))) r1s)
		  (ve3)))
	 (a2g (if (is> r2ss 0)
		  (* (- (/ mu2 (expt r2ss 3))) r2s)
		  (ve3))))
    (- (+ a1g a2g))))

(defun gradphi (res r12 mu1 mu2)
  (let ((omee (ve3 :c100 (/ (sqrt (+ mu1 mu2)) (sqrt (expt r12 3))))))
    (*x2 omee (*x2 omee res))))

(defun gradu (res r12 mu1 mu2)
  (+ (gradv res r12 mu1 mu2)
     (gradphi res r12 mu1 mu2)))

(defun sailcalcfideal (r1s gradu mu1 stype)
  (let ((r1suv (unitg r1s))
	(n (unitg gradu)))
    (cond ((is>= (dot r1s gradu) 0d0)
	   (let ((alpha (atan (/ (norme (*x2 r1suv gradu))
				 (dot r1suv gradu))))
		 (beta (* (/ (norme2 r1s) mu1)
			  (/ (dot gradu n)
			     (expt (dot r1suv n) stype)))))
	     (values beta alpha nil)))
	  (t (values 0 0 t)))))

(defun sailcalcideal (sailhash)
  (with-keys (r12 mu1 mu2 stype xmin xmax xsteps ymin ymax ysteps zmin zmax zsteps) sailhash
    (let* ((mu (/ mu2 (+ mu1 mu2)))
	   (re1 (ve3 :c1 (- (* r12 mu))))
	   (re2 (ve3 :c1 (* r12 (- 1 mu))))
	   (beta3d (make-array (list xsteps ysteps zsteps)))
	   (alpha3d (make-array (list xsteps ysteps zsteps))))
      (loop for ix below xsteps
	 for x = (+ xmin (gref re2 #b1) (* (/ (- xmax xmin) xsteps) ix))
	 do (loop for iy below ysteps
	       for y = (+ ymin (* (/ (- ymax ymin) ysteps) iy))
	       do (loop for iz below zsteps
		     for z = (+ zmin (* (/ (- zmax zmin) zsteps) iz))
		     for res = (ve3 :c1 x :c10 y :c100 z)
		     for gradu = (gradu res r12 mu1 mu2)
		     for r1s = (- res re1)
		     for (beta alpha err) = (multiple-value-list (sailcalcfideal r1s gradu mu1 stype))
		     do (setf (aref beta3d ix iy iz) beta)
		       (setf (aref alpha3d ix iy iz) alpha))))
      (values beta3d alpha3d))))

(defun sailcalcf (r1s gradu mu1 stype sigma* ref)
  (let* ((r1suv (unitg r1s))
	 (m (unitg gradu)))
    (cond ((is>= (dot r1s gradu) 0d0) ; gradu away from sun
	   (let* ((tanth (/ (norme (*x2 r1suv gradu))
			    (dot r1suv gradu)))
		  (temptanphi (* (/ (- 1 (expt ref 2))
				    (expt ref 2))
				 (expt tanth 2))))
	     (cond ((is> temptanphi 1) ; location invalid due to reflectivity limit on sail angle
		    (values 0 0 t))
		   ((is= tanth 0) ; r1s & m coincident: simpler equations
		    (let ((beta (* (/ (* 2 (norme2 r1s)) mu1)
				   (/ (norme gradu)
				      (+ 1 ref)))))
		      (values (/ sigma* beta) 0d0 nil)))
		   (t (let* ((tanphi (* (/ ref (* (+ 1 ref) tanth)) ; r1s & m not coincident: full equations
					(- 1 (sqrt (- 1 temptanphi)))))
			     (tanalpha (* (/ (+ 1 ref)
					     (- 1 ref))
					  tanphi))
			     (phi (atan tanphi))
			     (rmn (rotor (*o r1suv m) phi)) ; rotor to rotate m to n
			     (n (rot m rmn))
			     (beta (* (/ (* 2 (norme2 r1s))
					 mu1)
				      (/ (dot gradu n)
					 (* (+ 1 ref)
					    (expt (dot r1suv n) stype)))))
			     (sigma (/ sigma* beta)))
			(values sigma (atan tanalpha) nil))))))
	  (t (values 0 0 t)))))

(defun sailcalc (sailhash)
  (with-keys (r12 mu1 mu2 stype sigma* ref xmin xmax xsteps ymin ymax ysteps zmin zmax zsteps) sailhash
    (let* ((mu (/ mu2 (+ mu1 mu2)))
	   (re1 (ve3 :c1 (- (* r12 mu))))
	   (re2 (ve3 :c1 (* r12 (- 1 mu))))
	   (sigma3d (make-array (list xsteps ysteps zsteps)))
	   (alpha3d (make-array (list xsteps ysteps zsteps))))
      (loop for ix below xsteps
	 for x = (+ xmin (gref re2 #b1) (* (/ (- xmax xmin) xsteps) ix))
	 do (loop for iy below ysteps
	       for y = (+ ymin (* (/ (- ymax ymin) ysteps) iy))
	       do (loop for iz below zsteps
		     for z = (+ zmin (* (/ (- zmax zmin) zsteps) iz))
		     for res = (ve3 :c1 x :c10 y :c100 z)
		     for gradu = (gradu res r12 mu1 mu2)
		     for r1s = (- res re1)
		     for (sigma alpha err) = (multiple-value-list (sailcalcf r1s gradu mu1 stype sigma* ref))
		     do (setf (aref sigma3d ix iy iz) sigma)
		       (setf (aref alpha3d ix iy iz) alpha))))
      (values sigma3d alpha3d))))

(defun plot-sail-xz (sailhash sigma3d alpha3d levels)
  (with-keys (r12 mu1 mu2 stype ref rc pc xmin xmax xsteps ymin ymax ysteps zmin zmax zsteps) sailhash
    (let ((sigma2d (make-array (list xsteps zsteps)))
          (alpha2d (make-array (list xsteps zsteps)))
	  (alphalevels (map 'vector #'identity
			    (loop for d = 10 then (+ d 10)
			       while (< d 90)
			       collect (* d (/ pi 180))))))
      (dotimes (i xsteps)
        (dotimes (j zsteps)
          (setf (aref sigma2d i j) (aref sigma3d i 0 j))
          (setf (aref alpha2d i j) (aref alpha3d i 0 j))))
      (let ((sigma (new-contour-plot sigma2d :x-min xmin :x-max xmax :y-min zmin :y-max zmax :fill-type :none :contour-levels levels))
	    (alpha (new-contour-plot alpha2d :x-min xmin :x-max xmax :y-min zmin :y-max zmax :fill-type :none :contour-levels alphalevels))
            (w (basic-window)))
        (add-plot-to-window w sigma)
;;        (add-plot-to-window w alpha)
	(edit-window-axis w :x :properties (append (list :minor-tick-grid :major-tick-grid) (get-window-axis-properties w :x)))
	(edit-window-axis w :y :properties (append (list :minor-tick-grid :major-tick-grid) (get-window-axis-properties w :x)))
        (render w "xwin")))))

(defun export-gnuplot-xz (filename sailhash data3d &optional (iy 0))
  "Export 3D data (sigma or alpha) to Gnuplot compatible surface data file."
  (with-keys (r12 mu1 mu2 stype ref rc pc xmin xmax xsteps ymin ymax ysteps zmin zmax zsteps) sailhash
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (loop for ix below xsteps
	 for x = (+ xmin (* ix (- xmax xmin)))
	 do (loop for iz below zsteps
	       for z = (+ zmin (* iz (- zmax zmin)))
	       for dataxz = (aref data3d ix iy iz)
	       do (format s "~&~a~%" (substitute #\E #\d (format nil "~a ~a ~a" x z dataxz))))
	 do (format s "~%")))))

(defparameter *sun-earth-sail* ; units: MKS
  (make-hash :r12 1.49597870691d11 ; 1 AU
             :mu1 1.3271244d20 ; solar grav param
             :mu2 3.986d14 ; Earth grav param
             :stype 2 ; flat sail: 2, compound sail: 1
             :ref 0.85d0 ; reflectivity
	     :sigma* 1.53 ; g/m^2
             :rc 1.49597870691d11
             :pc 4.563d-6
             :xmin -3d9
             :xmax 1.6d9
             :xsteps 200
             :ymin 0d0
             :ymax 0d0
             :ysteps 1
             :zmin -2d9
             :zmax 2d9
             :zsteps 100))

(defparameter *solar-wind-nh-comm* ; solar wind with northern hemisphere comm
  (make-hash :r12 1.49597870691d11 ; 1 AU
             :mu1 1.3271244d20 ; solar grav param
             :mu2 3.986d14 ; Earth grav param
             :stype 2 ; flat sail: 2, compound sail: 1
             :ref 0.85d0 ; reflectivity
	     :sigma* 1.53 ; g/m^2
             :rc 1.49597870691d11
             :pc 4.563d-6
             :xmin -2d9
             :xmax -1.5d9
             :xsteps 100
             :ymin 0d0
             :ymax 0d0
             :ysteps 1
             :zmin 0.6d9
             :zmax 0.8d9
             :zsteps 100
	     :sigmalevels #(1d-6 26 44 46 48)))