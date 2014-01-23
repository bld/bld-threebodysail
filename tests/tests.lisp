(ql:quickload :bld-threebodysail)
(ql:quickload :fiveam)
(in-package :bld-threebodysail)
(use-package :fiveam)
(def-suite bld-threebodysail)
(in-suite bld-threebodysail)

;; Vector tests
(test isequal
  (is (isequal (* 0.9 *double-eps*) 0d0))
  (is (isequal (* -0.9 *double-eps*) 0d0))
  (is (not (isequal (* 1.1 *double-eps*) 0d0))))
(test isgt
  (is (isgt (* 1.1 *double-eps*) 0d0))
  (is (not (isgt (* 0.9 *double-eps*) 0d0))))
(test isgte
  (is (isgte (* 0.9 *double-eps*) 0d0))
  (is (isgte (* -0.9 *double-eps*) 0d0))
  (is (isgte (* 1.1 *double-eps*) 0d0))
  (is (not (isgte (* -1.1 *double-eps*) 0d0))))
(test +
  (is (every #'= (+ #(1 2 3) #(4 5 6)) #(5 7 9))))
(test -
  (is (every #'= (- #(1 2 3) #(4 5 6)) #(-3 -3 -3))))
(test *
  (is (every #'= (* #(1 2 3) 4) #(4 8 12)))
  (is (every #'= (* 4 #(1 2 3)) #(4 8 12))))
(test dot
  (is (= (dot #(1 2 3) #(1 2 3)) 14)))
(test norm2
  (is (isequal (norm2 #(1d0 2d0 3d0)) (sqrt 14d0))))
(test /
  (is (every #'= (/ #(1 2 3) 4) #(1/4 1/2 3/4))))
(test cross
  (is (every #'= (cross #(1 0 0) #(0 1 0)) #(0 0 1)))
  (is (every #'= (cross #(0 1 0) #(1 0 0)) #(0 0 -1)))
  (is (every #'= (cross #(0 1 0) #(0 0 1)) #(1 0 0))) 
  (is (every #'= (cross #(0 0 1) #(0 1 0)) #(-1 0 0)))
  (is (every #'= (cross #(0 0 1) #(1 0 0)) #(0 1 0)))
  (is (every #'= (cross #(1 0 0) #(0 0 1)) #(0 -1 0)))
  (is (every #'= (cross #(1 1 1) #(1 1 1)) #(0 0 0))))
(test coincident
  (is (coincident #(1 1 1) #(-2 -2 -2)))
  (is (not (coincident #(1 1 1) #(1 1 2)))))
(test codirectional
  (is (codirectional #(1 1 1) #(2 2 2)))
  (is (not (codirectional #(1 1 1) #(-1 -1 -1)))))

;; graduf
(test graduf
  (is (vectorp (graduf #(1 0 0) 1 1 0.01)))
  (is (= (length (graduf #(1 0 0) 1 1 0.01))) 3))
;; NEED A LOT MORE HERE

;; sailcalcf
;; NEED SOMETHING HERE

;; sailcalcloop
;; NEED SOMETHING HERE
