function LG = lgmode(p,l,r,phi,z,theta)
% LGMODE calculates LG mode amplitude at z = 0
%
% A = LGMODE(p,l,r,phi) calculates the LG mode amplitude for mode [p,l]
% at locations given in polar coordinates [r, phi].
% r is in units of the beam width; r and phi can be matrices of equal size.
%
% A = LGMODE(p,l,r,phi,z) computes the modes with z as well.
%
% A = LGMODE(p,l,r,phi,z,theta) scales the beam waist according to the
% beam convergence angle theta (in degrees): w0=1/tan(theta).
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:lgmode:move', ...
    'This file will move to ott.utils.lgmode');

w=paraxial_beam_waist(2*p+abs(l)); %Beam waist in normalized units.

if nargin<6
    theta=atan(w.^2/pi)*180/pi;
end
if nargin<5
	z=0;
end

k=2*pi;
w0=w/pi/tan(abs(theta/180*pi));
zr=pi*w0^2;
w_z=w0*sqrt(1+(z/zr).^2);
R=1./z*(z.^2+zr.^2);
R(z==0)=inf;
psi=atan2(z,zr);

extFac=exp(1i*(l*phi-k*z-r.^2/2./R));

r=r./w_z;

LG = sqrt(2*factorial(p)/(pi*factorial(p+abs(l)))) * (sqrt(2)*r).^abs(l) .* laguerre(p,abs(l),2*r.^2) ...
    .* extFac .* exp(-r.^2) .* exp(1i * (2*p + abs(l) + 1) * psi )./w_z;

return
