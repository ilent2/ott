function [nn,mm,a,b]=bsc_bessel(nmax,theta,lmode,Etheta,Ephi);
% bsc_bessel.m - calculate a bessel beam, and bessel-like beam with OAM.
%
% Usage:
% [n,m,a,b] = bsc_bessel(nmax,theta,[Ex,Ey]);
% OR
% [n,m,a,b] = bsc_bessel(nmax,theta,Etheta,Ephi);
% OR
% [n,m,a,b] = bsc_bessel(nmax,theta,lmode,[Ex,Ey]);
% OR
% [n,m,a,b] = bsc_bessel(nmax,theta,lmode,Etheta,Ephi);
%
% inputs are column vectors except for [Ex,Ey] which are jones vectors
% stacked in a column.
%
% lmode - applies orbital angualr momentum to the beam. 0 = ordinary
% bessel.
%
% This code uses bsc_plane to generate the basis modes.
%
% PACKAGE INFO

%% debug code
% 
% nmax=2;
% theta=pi-pi*linspace(.2,.4,1).';
% lmode=[1]%;1];
% 
% %Radial azimuthal:
% % Etheta=[1;1]/sqrt(2);
% % Ephi=[1i;1i]/sqrt(2);
% 
% %Uniform:
% Etheta=[1,1]%;1,-1i];
% %

if or(size(theta,2)~=1,size(lmode,1)~=size(theta,1))
    error('ott:bsc_bessel:badinput','Input is not valid. Inputs must be column vectors.')
end

if nargin==3
    %this is sufficient.
    Etheta=lmode;
    lmode=0*theta;
end

if nargin==4
    %two cases distinguishing only is ETheta is NOT wide.
    if size(Etheta,2)==1
        Ephi=Etheta;
        Etheta=lmode;
        lmode=0*theta;
    else
        %do nothing
    end
end

%nargin==5 is default behaviour

%% preamble
nTheta=length(theta);
indexes=[1:length(theta)].';

%% uniform modifications.
if size(Etheta,2)==2
    lmode=lmode+[-1,1]; %left is -1i for Ephi, right is +1i for Ephi;
    lmode=lmode(:);
    indexes=[[1:size(Etheta,1)].';[1:size(Etheta,1)].']; %packing for theta uniform
    
    %need to record polarisation and modify Etheta and Ephi.
    polarisation_weights=([1,1i;1,-1i]*(Etheta.')).'/2;
    Ephi=[-1i*ones(size(Etheta,1),1);1i*ones(size(Etheta,1),1)].*polarisation_weights(:);
    Etheta=[ones(size(Etheta,1),1);ones(size(Etheta,1),1)].*polarisation_weights(:).*sign(cos(theta(indexes)));
    
end

%% calculate the mode indices we are going to find.
% [nn,mm]=combined_index([1:nmax*(nmax+2)].');
a = zeros((nmax*(nmax+2)),nTheta);
b = zeros((nmax*(nmax+2)),nTheta);

for n = 1:nmax
    
    ci_index_m=find(abs(lmode)<=n);
    indt=n+lmode(ci_index_m)+1;
    
    %power normalisation.
    Nn = 1/sqrt(n*(n+1));
    
    %Generate the farfield components of the VSWFs
    [~,dtY,dpY] = spharm(n,theta);
         
    %slow indexing.
    szA=sub2ind(size(a),(n-1)*(n+1)+indt,indexes(ci_index_m));
    szY=sub2ind(size(dtY),indexes(ci_index_m),indt);
    
    dtY=dtY(:);
    dpY=dpY(:);
    
    %equivalent to dot((1i)^(n+1)*C,E);
    a(szA) = 4*pi*Nn*(-1i)^(n+1)*(conj(dpY(szY)).*Etheta(ci_index_m) - conj(dtY(szY)).*Ephi(ci_index_m));
    %equivalent to dot((1i)^(n)*B,E);
    b(szA) = 4*pi*Nn*(-1i)^(n)  *(conj(dtY(szY)).*Etheta(ci_index_m) + conj(dpY(szY)).*Ephi(ci_index_m));
    
end

%%

if nargout==4
    ci=find(any(abs(a)|abs(b),2));
    [nn,mm]=combined_index(ci);
    a=a(ci,:);
    b=b(ci,:);
end

if nargout==2
    nn=sparse(a);
    mm=sparse(b);
end

end
    