(TeX-add-style-hook
 "ps1"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10")
   (LaTeX-add-labels
    "eq:size"
    "eq:switch"
    "eq:xstar1")))

