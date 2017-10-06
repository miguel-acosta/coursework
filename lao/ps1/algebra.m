clear all
syms phix phip beta kx kp p theta gamma

PHIX = solve(phix == (1-beta + beta * phix)*(kx/(kx+kp)));

PHIP = solve(phip == (1-beta-beta*PHIX)*kp/(kx+kp) + beta * phip);

R    = (1-beta)*theta + beta*(PHIX*theta + PHIP*p);

deltap = (kx +kp)/(gamma*(1-beta+beta*PHIX)^2) - kp/(gamma*(1-beta+beta*PHIX));
deltax = -kx/(gamma*(1-beta + beta * PHIX));
