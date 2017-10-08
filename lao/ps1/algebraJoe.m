clear all
syms phix phip beta kx ke p theta gamma Phi1 nu Phi2

Phi1 = (1-beta)*kx/((1-beta)*kx + ke)

phi1 = kx*(1+(Phi1-1)*beta)/(kx+ke)

phi2 = ke*(1+(Phi1-1)*beta)/(kx+ke) + beta*Phi2

nu = ((1-beta)*ke + beta *Phi1 *ke)/(kx+ke)
Phi2 = nu/(1-beta)

