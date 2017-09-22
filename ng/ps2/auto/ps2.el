(TeX-add-style-hook
 "ps2"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "q1"
    "ogmm_ac"
    "k3_pretty"
    "k4_pretty"
    "ogmm_mom"
    "aux_ols_full"
    "aux_ols"
    "aux_ols_trans"
    "smm"
    "ii"
    "article"
    "art10")
   (LaTeX-add-labels
    "eq:g0g1"
    "eq:g2")))

