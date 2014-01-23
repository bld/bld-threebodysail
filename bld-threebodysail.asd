(asdf:defsystem :bld-threebodysail
  :author "Benjamin L. Diedrich <bldiedrich@gmail.com>"
  :license "MIT"
  :description "Calculates sun-planet-sail restricted 3-body equilibrium points, giving the sail performance and pointing required at points in the defined space"
  :depends-on (:bld-gen :bld-ga :bld-e3 :bld-utils)
  :components ((:module "src"
			:components
			((:file "package")
			 (:file "threebodysail" :depends-on ("package"))))))
