(in-package :bld-threebodysail)

(defparameter *double-eps* 1d-10)

(defun is= (a b)
  "Test if two numbers equal within tolerance"
  (< (abs (- (abs a) (abs b))) *double-eps*))

(defun is> (a b)
  "Test if A greater than B within tolerance"
  (> (- a b) *double-eps*))

(defun is>= (a b)
  "Test that A >= B within tolerance"
  (or (is> a b) (is= a b)))

(defmethod coincident ((a ve3) (b ve3))
  "Test if two vectors are coincident within tolerance"
  (is= (norme2 (*x2 a b)) 0d0))

(defmethod codirectional ((a ve3) (b ve3))
  "Test if two vectors point in the same direction within tolerance"
  (and (coincident a b)
       (is> (scalar (*i a b)) 0d0)))

(defmethod dot ((a ve3) (b ve3))
 "Dot product of two vectors"
  (scalar (*i a b)))

(defun gradv (res r12 mu1 mu2)
  "Gradient of gravitational potential"
  (let* ((mu (/ mu2 (+ mu1 mu2)))
	 (re1 (ve3 :e1 (* -1 r12 mu)))
	 (re2 (ve3 :e1 (* r12 (- 1 mu))))
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
 "Gradient of centripetal potential"
  (let ((omee (ve3 :e3 (/ (sqrt (+ mu1 mu2)) (sqrt (expt r12 3))))))
    (*x2 omee (*x2 omee res))))

(defun gradu (res r12 mu1 mu2)
  "Gradient of combined gravitational/centripetal potential"
  (+ (gradv res r12 mu1 mu2)
     (gradphi res r12 mu1 mu2)))

(defun sailcalcfideal (r1s gradu mu1 stype)
  "Calculate loading & cone angle for ideal sail at sun-sail vector"
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
  "Calculate 3D arrays of loading & cone angle over space & parameters defined in hash table"
  (with-keys (r12 mu1 mu2 stype xmin xmax xsteps ymin ymax ysteps zmin zmax zsteps) sailhash
    (let* ((mu (/ mu2 (+ mu1 mu2)))
	   (re1 (ve3 :e1 (- (* r12 mu))))
	   (re2 (ve3 :e1 (* r12 (- 1 mu))))
	   (beta3d (make-array (list xsteps ysteps zsteps)))
	   (alpha3d (make-array (list xsteps ysteps zsteps))))
      (loop for ix below xsteps
	 for x = (+ xmin (gref re2 :e1) (* (/ (- xmax xmin) xsteps) ix))
	 do (loop for iy below ysteps
	       for y = (+ ymin (* (/ (- ymax ymin) ysteps) iy))
	       do (loop for iz below zsteps
		     for z = (+ zmin (* (/ (- zmax zmin) zsteps) iz))
		     for res = (ve3 :e1 x :e2 y :e3 z)
		     for gradu = (gradu res r12 mu1 mu2)
		     for r1s = (- res re1)
		     for (beta alpha err) = (multiple-value-list (sailcalcfideal r1s gradu mu1 stype))
		     do (setf (aref beta3d ix iy iz) beta)
		       (setf (aref alpha3d ix iy iz) alpha))))
      (values beta3d alpha3d))))

(defun sailcalcf (r1s gradu mu1 stype sigma* ref)
  "Calculate sail loading & cone angle at vector from 1st body (sun) to sail, gravitational parameter of sun, sail stype, critical sail loading, and reflection."
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
			     (n (rotateg m rmn))
			     (beta (* (/ (* 2 (norme2 r1s))
					 mu1)
				      (/ (dot gradu n)
					 (* (+ 1 ref)
					    (expt (dot r1suv n) stype)))))
			     (sigma (/ sigma* beta)))
			(values sigma (atan tanalpha) nil))))))
	  (t (values 0 0 t)))))

(defun sailcalc (sailhash)
  "Generate 3D arrays of sail loading & cone angles required for equilibrium given hash table of sail problem parameters"
  (with-keys (r12 mu1 mu2 stype sigma* ref xmin xmax xsteps ymin ymax ysteps zmin zmax zsteps pc) sailhash
    (let* ((mu (/ mu2 (+ mu1 mu2)))
	   (re1 (ve3 :e1 (- (* r12 mu))))
	   (re2 (ve3 :e1 (* r12 (- 1 mu))))
	   (sigma3d (make-array (list (1+ xsteps) (1+ ysteps) (1+ zsteps))))
	   (alpha3d (make-array (list (1+ xsteps) (1+ ysteps) (1+ zsteps))))
	   (ac3d (make-array (list (1+ xsteps) (1+ ysteps) (1+ zsteps)))))
      (loop for ix upto xsteps
	 for x = (+ xmin (gref re2 :e1) (* (/ (- xmax xmin) xsteps) ix))
	 do (loop for iy upto ysteps
	       for y = (+ ymin (* (/ (- ymax ymin) ysteps) iy))
	       do (loop for iz upto zsteps
		     for z = (+ zmin (* (/ (- zmax zmin) zsteps) iz))
		     for res = (ve3 :e1 x :e2 y :e3 z)
		     for gradu = (gradu res r12 mu1 mu2)
		     for r1s = (- res re1)
		     for (sigma alpha err) = (multiple-value-list (sailcalcf r1s gradu mu1 stype sigma* ref))
		     for ac = (if (> sigma 0)
				  (* (/ (* 2 ref pc) (/ sigma 1000)) 1000)
				  1e6)
		     do (setf (aref sigma3d ix iy iz) sigma
			      (aref alpha3d ix iy iz) alpha
			      (aref ac3d ix iy iz) ac))))
      (values sigma3d alpha3d ac3d))))

(defun export-gnuplot-xz (filename sailhash data3d &key (iy 0) (xyzscale 1d-9) (datascale 1d0))
  "Export 3D data (sigma or alpha) to Gnuplot compatible surface data file."
  (with-keys (r12 mu1 mu2 stype ref rc pc xmin xmax xsteps ymin ymax ysteps zmin zmax zsteps) sailhash
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (loop for ix upto xsteps
	 for x = (* xyzscale (+ xmin (* (/ (- xmax xmin) xsteps) ix)))
	 do (loop for iz upto zsteps
	       for z = (* xyzscale (+ zmin (* (/ (- zmax zmin) zsteps) iz)))
	       for dataxz = (* datascale (aref data3d ix iy iz))
	       do (format s "~f ~f ~f~%" x z dataxz))
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

(defparameter *sun-earth-sail-subl1* ; units: MKS
  (make-hash :r12 1.49597870691d11 ; 1 AU
             :mu1 1.3271244d20 ; solar grav param
             :mu2 3.986d14 ; Earth grav param
             :stype 2 ; flat sail: 2, compound sail: 1
             :ref 0.85d0 ; reflectivity
	     :sigma* 1.53 ; g/m^2
             :rc 1.49597870691d11
             :pc 4.563d-6
             :xmin -3d9
             :xmax 0d9
             :xsteps 200
             :ymin 0d0
             :ymax 0d0
             :ysteps 1
             :zmin -2d9
             :zmax 2d9
             :zsteps 100))
