function [rtp,n,ds]=axisym_boundarypoints(Nmax,rho,z);
% axisym_boundarypoints.m : Calculates boundary points for surface integral
%                           given a perimeter in r, theta_sp and the
%                           normals n. Is in spherical coordinates.
%
% Usage:
% [rhoout,zout]=axisym_boundarypoints(rho,z,n);
%
% PACKAGE INFO

% %%%%test values
% ntheta=500;
% rho=[0,1,1,0];
% z=[1,1,-1,-1];
% % rho=sin(linspace(0,pi,100));
% % z=cos(linspace(0,pi,100));
% %%%%

%column vectors only, thank you.
rho=rho(:);
z=z(:);

ntheta=(Nmax+2)*1024;

%add up line segments to get distance
s=sqrt((rho(2:end)-rho(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2);

ds=sum(s)/ntheta;

%Don't ask me how this works. It does. It's simple algebra in the end...
%and yes, it can be done more intelligently.
zout=zeros(ntheta,1);
rhoout=zout;
nxyz=zeros(ntheta,3);

sdeficit=0;
ncum=0;

for ii=2:length(rho);
    N=s(ii-1)/ds;
    Nused=round(N+sdeficit);
    
    nc=[-(z(ii)-z(ii-1)),0,(rho(ii)-rho(ii-1))];
    nc=nc/norm(nc,2);
    
    if Nused>=1
        drho=(rho(ii)-rho(ii-1))/N*ones(Nused,1);
        
        rhot=cumsum(drho)-drho/2-sdeficit*drho(1);
        rhoout([ncum+[1:Nused]])=rho(ii-1)+rhot;
        
        dz=(z(ii)-z(ii-1))/N*ones(Nused,1);
        
        zt=cumsum(dz)-dz/2-sdeficit*dz(1);
        zout([ncum+[1:Nused]])=z(ii-1)+zt;
        
        nxyz([ncum+[1:Nused]],:)=repmat(nc,[length(zt),1]);
        
        sdeficit=(N-Nused+sdeficit);
    else
        sdeficit=sdeficit+N;
    end
    
    ncum=ncum+Nused;
    
end

%converts the cylindrical coordinates into spherical coordinates
[n,rtp]=xyzv2rtpv(nxyz,[rhoout,zeros(size(rhoout)),zout]);

dst=zeros(length(rhoout),3);

%calcultes area elements
dst(2:end-1,1)=(rhoout(3:end)-rhoout(1:end-2))/2;
dst(2:end-1,3)=(zout(3:end)-zout(1:end-2))/2;
dst(1,1)=(mean(rhoout(1:2))-rho(1));
dst(1,3)=(mean(zout(1:2))-z(1));
dst(end,1)=(rho(end)-mean(rhoout(end-1:end)));
dst(end,3)=(z(end)-mean(zout(end-1:end)));

%a general axisymmetric conic region has the following area (sans factor of 2*pi):
ds=(rtp(:,1).*sqrt(sum(abs(dst).^2,2)).*sin(rtp(:,2))); 

