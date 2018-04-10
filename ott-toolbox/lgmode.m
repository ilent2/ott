function LG = lgmode(p,l,r,phi,z,theta)
% LGMODE calculates LG mode amplitude at z = 0
%
% A = LGMODE(p,l,r,phi) calculates the LG mode amplitude for mode [p,l]
% at locations given in polar coordinates [r, phi].
% r is in units of the beam width; r and phi can be matrices of equal size.
%
% A = LGMODE(p,l,r,phi,theta) scales the beam waist according to the
% beam convergence angle theta (in degrees): w0=1/tan(theta).
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin<6
    theta=atan(1/pi)*180/pi;
end
if nargin<5
	z=0;
end

w = 1.; %Beam waist in normalized units.

% if (l ~= 0)
%     invL=1./abs(l );
%     zz = exp(-(abs(l )+2.)*invL);
%     w=-(1.+2*sqrt(invL)+invL); %This is a really good starting guess. It converges within 3 iterations for l=1:10000+
%     
%     w0=-w;
%     
%     while (abs(w-w0)>0.00001)
%         w0=w;
%         expw = exp(w);
%         
%         w=w0-(w0*expw+zz)/(expw+w0*expw); %Newton's rule... Usually this kind of method would find the real root i.e. W_0(z)... This finds W_{-1}(z) local to the beam waist of an LG beam.
%         
%     end
%     
%     w = sqrt(-abs(l )/2.*w); %Beam waist in normalized units
%     
% end

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
