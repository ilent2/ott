function [fx,fy,fz,tx,ty,tz]=force_torque_farsund(n,m,a,b,p,q)
% force_torque_farsun.m : Calculates the force and torque in a 3D orthogonal
%                         space. If the beam shape coefficients are in the
%                         original coordinates, this outputs the force and
%                         torque in 3D carteisan coordinates.
%
% USAGE:
%
% [fx,fy,fz,tx,ty,tz] = force_torque_farsund(n,m,a,b,p,q)
%
% OR
% 
% [fxyz,txyz] = force_torque_farsund(n,m,a,b,p,q)
%
% OR
%
% [fxyz] = force_torque_farsund(n,m,a,b,p,q)
%
% where n,m are the degree and order of modes contained in the system. They
% can be the truncated n and m's.
%
% a,b,p,q are the incident modes and scattered modes respectively. These 
% must be full or sparse of the appropriate size!!!!!
%
% This uses mathematical result of Farsund et al., 1996, in the form of
% Chricton and Marsden, 2000, and our standard T-matrix notation S.T. 
% E_{inc}=sum_{nm}(aM+bN);
%
% PACKAGE INFO

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
fy=-imag(fxy);

tz=sum(m.*(a.*conj(a)+b.*conj(b)-p.*conj(p)-q.*conj(q))); %this has the correct sign... farsund. modes match.

txy=sum(sqrt((n-m).*(n+m+1)).*(a.*conj(amp1)+b.*conj(bmp1)-p.*conj(pmp1)-q.*conj(qmp1))); %this has the correct sign... farsund. modes match.
tx=real(txy);
ty=-imag(txy);

if nargout <= 2
    fx=[fx;fy;fz];
    fy=[tx;ty;tz];
end