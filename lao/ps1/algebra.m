clear all
syms phix phip beta kx kp p theta gamma x lambda ke

PHIX = solve(phix == (1-beta + beta * phix)*(kx/(kx+kp)));

PHIP = solve(phip == (1-beta+beta*PHIX)*kp/(kx+kp) + beta * phip);

R    = (1-beta)*theta + beta*(PHIX*theta + PHIP*p);

ER = (1-beta+beta*PHIX)*(kx*x + kp*p)/(kx+kp) + beta *PHIP*p;

VR = (1-beta + beta*PHIX)^2/(kx+kp);

kj = (p-ER)/(gamma*VR);

m = lambda - gamma*(ke/(lambda^2) + kx)*(1-beta)/(kx*ke/(lambda^2) + (1-beta)*kx);

pretty(simplify(-diff(m,ke)/diff(m,lambda)))
pretty(simplify(-diff(m,gamma)/diff(m,lambda)))

