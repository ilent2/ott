function [E,H] = farfield(n,m,a,b,p,q,theta,phi)
%FARFIELD finds far field given VSWF expansion coefficients
%
% [E,H] = FARFIELD(n,m,a,b,p,q,theta,phi) calculates far fields.
% Each row of E,H is the field (in spherical coordinates) in the
% (theta,phi) direction (assuming a distance scaling factor of kr)
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('internal');

[theta,phi] = matchsize(theta,phi);

[theta_new,~,indY]=unique(theta);
[phi_new,~,indP]=unique(phi);

E = zeros(length(theta),3);
Etheta=zeros(length(theta),1);
Ephi=zeros(length(theta),1);

H = zeros(length(theta),3);
Htheta=zeros(length(theta),1);
Hphi=zeros(length(theta),1);

a = threewide(a);
b = threewide(b);
p = threewide(p);
q = threewide(q);

for nn = 1:max(n)
    
    vv=find(n==nn);
    [Y,Ytheta,Yphi] = spharm(nn,m(vv),theta_new,zeros(size(theta_new)));
    
    [M,PHI]=meshgrid(m(vv),phi_new);
    
    expimphi=repmat(exp(1i*M.*PHI),[1,3]);
    
    %this makes the vectors go down in m for n. has no effect if old version
    %code.
    Nn = 1/sqrt(nn*(nn+1));
    
    for ii=1:length(vv)
        index=nn*(nn+1)+m(vv(ii));
        
        TEMP=Nn*(((1i)^(nn+1) * a(index)+(-1i)^(nn+1) * p(index)).*Yphi(:,ii)+((1i)^nn * b(index)+(-1i)^nn * q(index)).*Ytheta(:,ii));
        Etheta=Etheta+TEMP(indY).*expimphi(indP,ii);
        
        TEMP=Nn*(-((1i)^(nn+1) * a(index)+(-1i)^(nn+1) * p(index)).*Ytheta(:,ii)+((1i)^nn * b(index)+(-1i)^nn * q(index)).*Yphi(:,ii));
        Ephi=Ephi+TEMP(indY).*expimphi(indP,ii);
        
        if nargout>1
            TEMP=Nn*(((1i)^(nn+1) * b(index)+(-1i)^(nn+1) * q(index)).*Yphi(:,ii)+((1i)^nn * a(index)+(-1i)^nn * p(index)).*Ytheta(:,ii));
            Htheta=Htheta+TEMP(indY).*expimphi(indP,ii);
            
            TEMP=Nn*(-((1i)^(nn+1) * b(index)+(-1i)^(nn+1) * q(index)).*Ytheta(:,ii)+((1i)^nn * a(index)+(-1i)^nn * p(index)).*Yphi(:,ii));
            Hphi=Hphi+TEMP(indY).*expimphi(indP,ii);
        end
        
    end
    
    
end
E=[zeros(size(Etheta)),Etheta,Ephi];
H=[zeros(size(Htheta)),Htheta,Hphi];

% SI-ify units of H
H = H * -1i;

ott_warning('external');

return
