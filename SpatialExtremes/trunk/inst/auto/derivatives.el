(TeX-add-style-hook "derivatives"
 (lambda ()
    (LaTeX-add-labels
     "sec:smith-char"
     "eq:smith"
     "sec:usefull-quantities"
     "eq:1"
     "sec:density-computation"
     "eq:2"
     "eq:6"
     "eq:3"
     "eq:4"
     "eq:5"
     "eq:10"
     "eq:7"
     "eq:9"
     "eq:11"
     "sec:gradient-computation"
     "sec:schlather-char")
    (TeX-run-style-hooks
     "amsmath"
     "a4wide"
     "latex2e"
     "art10"
     "article")))

