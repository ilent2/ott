function HG=hgmode(m,n,x,y,z,theta);
% HGMODE calcualtes HG mode amplitude at z = 0
%
% A = HGMODE(m,n,x,y) calcualtes HG mode amplitude for mode [m,n] at
% location [x,y].
%
% A = HGMODE(m,n,x,y,z) calcualtes HG mode amplitude for mode [m,n] at
% location [x,y,z].
%
% A = HGMODE(m,n,x,y,z,theta) scales the beam waist according to the
% beam convergence angle theta (in degrees): w0=1/tan(theta).

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

w=ott.utils.paraxial_beam_waist(m+n); %Beam waist in normalized units.

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

r=sqrt(x.^2+y.^2);

extFac=exp(-1i*(k*z+r.^2/2./R));

r=r./w_z;
x=x./w_z;
y=y./w_z;

HG = exp(1i*(n + m + 1)*psi)*(2/pi)^(1/2).*...
    (1./(2^(n+m)*factorial(n)*factorial(m)*w_z)).^(1/2).*...
    ott.utils.hermite(m, (sqrt(2)*x)).*ott.utils.hermite(n, (sqrt(2)*y)).* ...
    exp( - r.^2).*extFac./sqrt(w_z);
