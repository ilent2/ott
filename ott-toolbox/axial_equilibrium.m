function [z,kz] = axial_equilibrium(T,a,b,z)
%AXIAL_EQUILIBRIUM find equilibrium position and stiffness along beam axis
%
% [z,kz] = AXIAL_EQUILIBRIUM(T,a,b) attempts to locate the equilibrium
% position for the T-matrix T in beam [a, b] starting with an initial
% guess at z = 0.
%
% [z,kz] = axial_equilibrium(..., initial_guess) specifies an initial guess.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('internal');

if nargin < 4
    z = 0;
end

szT=size(T)/2;
sza=length(a);

if max(szT)>sza %zero pad a,b
    anew=sparse(szT(1),1);
    bnew=anew;
    anew(1:sza)=a;
    bnew(1:sza)=b;
    a=anew;
    b=bnew;
    clear anew bnew;
    nmax=floor(sqrt(szT(1)));
else
    Tnew=sparse(2*szT(2),2*length(a));
    Tnew(:,1:szT(1))=T(:,1:szT(1));
    Tnew(:,sza+1:sza+szT(1))=T(:,szT(1)+1:end);
    nmax=floor(sqrt(sza));
    T=Tnew;
    clear Tnew;
end

equiv_ka = nmax2ka(nmax);

power=sqrt(sum(abs(a).^2+abs(b).^2));

a=a/power;
b=b/power;

[n,m]=combined_index(1:nmax^2+2*nmax);

%start with three points each side over a 1/8 wavelength:
zd=.25;
dz=zd/5;
zs0=[z-zd/2:dz:z+zd/2];
zs=zs0;

fz=zeros(size(zs));

chkflg=false;
jj=0;

p=zeros(length(a),1);
q=p;

while ~chkflg
    jj=jj+1;
    for ii=1:length(zs)
        if fz(ii)==0
            [A,B] = translate_z(nmax,zs(ii));
            a1 = ( A*a + B*b );
            b1 = ( A*b + B*a );
            
            pq=T*[a1;b1];
            
            p(1:length(pq)/2)=pq(1:length(pq)/2);
            q(1:length(pq)/2)=pq(length(pq)/2+1:end);
            
            %power has already been normalized to 1.
            [~,~,fz(ii),~,~,~]=forcetorque(n(:),m(:),a1(:),b1(:),p(:),q(:));
        end
    end
    
    dzf=diff(fz)./diff(zs);
    
    indf=find(fz<=0,1);
    indzf=find(dzf<=0,1);
    
    if length(indf)
        if indf==1
            zs=[zs0(1:end-1)-zd*jj,zs];
            fz=[zs0(1:end-1)*0,fz];
        elseif indf==length(fz)
            zs=[zs,zs0(2:end)+zd*jj];
            fz=[fz,zs0(2:end)*0];
        else
            chkflg=true;
        end
    else
        if indzf>0
            zs=[zs,zs0(2:end)+zd*jj];
            fz=[fz,zs0(2:end)*0];
        elseif indzf<0
            zs=[zs0(1:end-1)-zd*jj,zs];
            fz=[zs0(1:end-1)*0,fz];
        else
            zs=[zs0(1:end-1)-zd*jj,zs,zs0(2:end)+zd*jj];
            fz=[zs0(1:end-1)*0,fz,zs0(2:end)*0];
        end
    end
    if jj>5
        ott_warning('external');
        error('No stable equilibrium near z!')
    end
end
pl=polyfit(zs(max([1,indf-3]):min([length(zs),indf+3])),fz(max([1,indf-3]):min([length(fz),indf+3])),3);

% %testing
% plot(zs,fz,'.')
% hold on
% plot(zs,polyval(pl,zs),'r')
% hold off

rts=roots(pl);
rts=rts(imag(rts)==0);

dzf=polyval([3*pl(1),2*pl(2),pl(3)],rts);

rtsi=find(dzf<0);

if length(rtsi)
    [d,indz]=min(abs(rts(rtsi)-z));
    z=rts(rtsi(indz));
    kz=dzf(rtsi(indz));
else
    ott_warning('external');
    error('No stable equilibrium near z!')
end

ott_warning('external');

return
