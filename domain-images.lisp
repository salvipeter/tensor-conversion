(in-package :cl-nurbs-tests)

(defmethod compute-distance ((type (eql 'tensor)) points segments p dir)
  (let ((p0 (elt segments 1))
	(p1 (elt segments 2)))
    (if (eq dir 'd)
	(/ (point-line-distance p (list p0 p1))
	   (point-line-distance (first segments) (list p0 p1)))
	(let ((di-1 (point-line-distance p (list (elt segments 0) p0)))
	      (di+1 (point-line-distance p (list p1 (elt segments 3)))))
	  (/ di-1 (+ di-1 di+1))))))

#+nil
(let ((points (points-from-angles (uniform-angles 5))))
  (sliced-distance-function-test points '(s s nil nil nil) "/tmp/s-params.ps"
                                 :resolution 30 :density 0.1
                                 :distance-type 'tensor :color t))

#+nil
(let ((points (points-from-angles (uniform-angles 5))))
  (sliced-distance-function-test points '(d d nil nil nil) "/tmp/h-params.ps"
                                 :resolution 30 :density 0.1
                                 :distance-type 'tensor :color t))

;;; + rotate 198 degrees CW
