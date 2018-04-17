function [fx,fy,fz,tx,ty,tz,sx,sy,sz]=forcetorque(n,m,a,b,p,q)
% FORCETORQUE calculate force/torque/spin in a 3D orthogonal space
% If the beam shape coefficients are in the original coordinates,
% this outputs the force, torque and spin in 3D carteisan coordinates.
%
% [fxyz,txyz,sxyz] = FORCETORQUE(n,m,a,b,p,q) calculates the force,
% torque and spin using the incident a,b beam vectors and the scattered
% p,q beam vectors.
%
% Output is stored in [3, 1] column vectors.  If torque or spin are
% omitted, only force or force/torque are calculated.
%
% n,m are the degree and order of modes contained in the system. They
% can be the truncated n and m's.
% a,b,p,q are the incident modes and scattered modes respectively. These
% must be full or sparse of the appropriate size!!!!!
%
% [...] = FORCETORQUE(n,m,ab,pq) as above but with the incident/scattered
% beam vectors specified by single column vectors ab = [a;b] and pq = [p;q].
%
% [fx,fy,fz,tx,ty,tz,sx,sy,sz] = FORCETORQUE(...) unpacks the
% force/torque/spin into separate output arguments.
%
% This uses mathematical result of Farsund et al., 1996, in the form of
% Chricton and Marsden, 2000, and our standard T-matrix notation S.T.
% E_{inc}=sum_{nm}(aM+bN);
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin == 4
  % Split ab and pq vectors
  ab = a;
  pq = b;
  a = ab(1:length(ab)/2);
  b = ab(length(ab)/2+1:end);
  p = pq(1:length(pq)/2);
  q = pq(length(pq)/2+1:end);
elseif nargin == 6
  % Nothing to do
else
  error('ott:forcetorque:nargin', 'Invalid number of arguments');
end

ott_warning('internal');

fx=0;
fy=0;
fz=0;
tx=0;
ty=0;
tz=0;
sx=0;
sy=0;
sz=0;

nmax=max(n);

b=1i*b;
q=1i*q;

if 1
    p=a+2*p;
    q=b+2*q;
end

addv=zeros(2*nmax+3,1);

at=[a;addv];
bt=[b;addv];
pt=[p;addv];
qt=[q;addv];

ci=combined_index(n,m);

%these preserve order and number of entries!
np1=2*n+2;
cinp1=ci+np1;
cinp1mp1=ci+np1+1;
cinp1mm1=ci+np1-1;
cimp1=ci+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is for m+1... if m+1>n then we'll ignore!
kimp=(m>n-1);

anp1=at(cinp1);
bnp1=bt(cinp1);
pnp1=pt(cinp1);
qnp1=qt(cinp1);

anp1mp1=at(cinp1mp1);
bnp1mp1=bt(cinp1mp1);
pnp1mp1=pt(cinp1mp1);
qnp1mp1=qt(cinp1mp1);

anp1mm1=at(cinp1mm1);
bnp1mm1=bt(cinp1mm1);
pnp1mm1=pt(cinp1mm1);
qnp1mm1=qt(cinp1mm1);

amp1=at(cimp1);
bmp1=bt(cimp1);
pmp1=pt(cimp1);
qmp1=qt(cimp1);

amp1(kimp)=0;
bmp1(kimp)=0;
pmp1(kimp)=0;
qmp1(kimp)=0;

a=a(ci);
b=b(ci);
p=p(ci);
q=q(ci);

Az=m./n./(n+1).*imag(-(a).*conj(b)+conj(q).*(p)); %this has the correct sign... farsund. modes match.
Bz=1./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ... %.*n
    .*imag(anp1.*conj(a)+bnp1.*conj(b)-(pnp1).*conj(p) ...
    -(qnp1).*conj(q)); %this has the correct sign... farsund. modes match.

fz=2*sum(Az+Bz);

Axy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).*(conj(pmp1).*q - conj(amp1).*b - conj(qmp1).*p + conj(bmp1).*a); %this has the correct sign... farsund. modes match.
Bxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ... %sqrt(n.*)
    ( sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(pnp1mp1) + q.* conj(qnp1mp1) -a.*conj(anp1mp1) -b.*conj(bnp1mp1)) + ... %this has the correct sign... farsund. modes match.
    sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(p) + qnp1mm1.*conj(q) - anp1mm1.*conj(a) - bnp1mm1.*conj(b)) ); %this has the correct sign... farsund. modes match.

fxy=sum(Axy+Bxy);
fx=real(fxy);
fy=imag(fxy);

if nargout > 1
    tz=sum(m.*(a.*conj(a)+b.*conj(b)-p.*conj(p)-q.*conj(q))); %this has the correct sign... farsund. modes match.
    
    txy=sum(sqrt((n-m).*(n+m+1)).*(a.*conj(amp1)+b.*conj(bmp1)-p.*conj(pmp1)-q.*conj(qmp1))); %this has the correct sign... farsund. modes match.
    tx=real(txy);
    ty=imag(txy);
    
    if nargout > 2
        Cz=m./n./(n+1).*(-(a).*conj(a)+conj(q).*(q)-(b).*conj(b)+conj(p).*(p));
        Dz=-2./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ...
            .*real(anp1.*conj(b)-bnp1.*conj(a)-(pnp1).*conj(q) ...
            +(qnp1).*conj(p));
        
        sz = sum(Cz+Dz);
        
        Cxy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).*(conj(pmp1).*p - conj(amp1).*a + conj(qmp1).*q - conj(bmp1).*b);
        Dxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ...
            ( (sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(qnp1mp1) - q.* conj(pnp1mp1) -a.*conj(bnp1mp1) +b.*conj(anp1mp1))) + ...
            (sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(q) - qnp1mm1.*conj(p) - anp1mm1.*conj(b) + bnp1mm1.*conj(a))) );
        
        sxy=sum(Cxy+Dxy);
        sy=real(sxy);
        sx=imag(sxy);
    end
end

if nargout <= 3
    fx=[fx;fy;fz];
    fy=[tx;ty;tz];
    fz=[sx;sy;sz];
end

