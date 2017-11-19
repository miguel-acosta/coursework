syms muI s l Pn ci ch

Pd = 1-Pn;
muH = 1-muI;

eq1 = muI*exp(-s/l)/(Pn*exp(-s/l) + Pd * exp(-ci/l)) + ...
      muH/(Pn + Pd * exp(-ch/l)) ;

eq2 = muI*exp(-ci/l)/(Pn*exp(-s/l) + Pd*exp(-ci/l)) + ...
      muH*exp(-ch/l)/(Pn + Pd*exp(-ch/l)) ;

%simplify([eq1; eq2])
