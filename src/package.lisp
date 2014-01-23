(defpackage :bld-threebodysail
  (:use :cl :bld-e3 :bld-ga)
  (:shadowing-import-from bld-gen
			  + - * / expt
			  sin cos tan
			  atan asin acos
			  sinh cosh tanh
			  asinh acosh atanh
			  log exp sqrt abs
			  min max signum)
  (:import-from :bld-utils make-hash with-keys print-object)
  (:export graduf sailcalcf sailcalcloop sailcalc))
