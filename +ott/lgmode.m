function LG = lgmode(p,l,r,phi)
% lgmode -- LG mode at z = 0
%
% Usage:
% A = lgmode(p,l,r,phi);
% where
% r is in units of the beam width
% r and phi can be matrices of equal size
%
% PACKAGE INFO

LG = sqrt(2*factorial(p)/(pi*factorial(p+abs(l)))) * (sqrt(2)*r).^abs(l) .* laguerre(p,abs(l),2*r.^2) ...
    .* exp(i*l*phi) .* exp(-r.^2) * exp(i * (2*p + abs(l) + 1) * pi/2);

return
