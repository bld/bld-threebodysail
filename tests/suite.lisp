(in-package :bld-threebodysail)

(use-package :fiveam)

(def-suite bld-threebodysail)

#|
(defmethod asdf:perform ((op asdf:test-op) (system (eql (asdf:find-system :bld-threebodysail-tests))))
  (format t "~2&********************~@
                ** Starting tests **~@
                ********************~%")
  (run! 'bld-threebodysail)
  (format t "~2&*********************************************~@
                **            Tests finished               **~@
                *********************************************~@
                ** If there were any failures on your      **~@
                ** platform, please report them to me:     **~@
                **    (bldiedrich at gmail dot com)        **~@
                ** or just file a bugreport on github:     **~@
                ** github.com/bld/bld-threebodysail/issues **~@
                *****************************************~%"))
|#
