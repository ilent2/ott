function [nn,mm,a,b]=paraxial_to_bsc(NA,E_ff,polarisation,varargin);
% PARAXIAL_TO_BSC Takes complex amplitudes in a plane (rho \propto \theta)
% and maps it onto an angular grid with extent dermined by NA, and solves 
% the BSCs.
% NOTE: Default behaviour is to have the BFP radius \propto sin(theta).
%
% [n,m,a,b]=paraxial_to_bsc(NA,E_ff,polarisation,(optional));
%
% [a,b]=paraxial_to_bsc(NA,E_ff,polarisation,(optional));
%
% inputs:
% NA : numerical aperture is the extent of the complex amplitudes.
% E_ff : complex amplitudes of the diffracted pattern. NOTE: rectangular 
%        grids use the SHORTEST dimension. The "pixels" are assumed square.
% polarisation : standard jones vector for polarisation.
%
% (optional) includes:
% polarisation --- { 'azimuthal' | 'radial' }. default: off
% BFP mapping --- { {'sintheta'} | 'tantheta' | 'theta' }. 'sintheta' uses 
%  BFP \rho \prop sin(\theta). 'tantheta' has the back  aperture \rho \prop 
%  tan(\theta). 'theta' uses a r \prop theta convention.
% nmax --- {{ 'nmax',nmax }} if a cell is input with 'nmax' as the first
%  element, the second argument is nmax, by default: nmax=30.
% fitting grid --- {{'grid',ntheta,nphi}} cell input of number of theta
%  points, ntheta, and phi points, nphi. 
%  default : ntheta=2*(nmax+1),  nphi=2*(nmax+1).
% refractive index --- {{'ri',nMeidum}}, cell input of refractive index.
%  default : nMedium=1.33;
%
% outputs:
% n,m : mode indices
% a,b : BSCs for VSWF.
%
% NOTE: This current version will best work for "clean" beam modes, that
% is, the desired field AFTER spatial filtering (If modelling an SLM/DMD).
% In addition, a standard G&L algorithm will produce abberations in the 
% output BSC.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

verbose=0;
nMedium=1.33;
Nmax=30;
function_theta=2;
ra=false;
az=false;
if nargin==2
    polarisation=[1,0];
end

if nargin>3
    for ii=1:length(varargin)
        switch class(varargin{ii})
            case {'char'}
                switch lower(varargin{ii})
                    case 'azimuthal'
                        az=true;
                        polarisation=[0,1];
                    case 'radial'
                        ra=true;     
                        polarisation=[1,0];
                    case 'theta'
                        function_theta=1;
                    case 'sintheta'
                        function_theta=2;
                    case 'tantheta'
                        function_theta=0;
                    otherwise
                        warning('ott:SLM_image_to_bsc:input',['The input ''' varargin{ii} ''' is not recognised! The code probably isn''t doing what you want it to as a result!'])
                end
                
            case {'cell'}
                switch lower(varargin{ii}{1})
                    case 'nmax'
                        Nmax=varargin{ii}{2};
                    case 'grid'
                        ntheta=varargin{ii}{2};
                        nphi=varargin{ii}{3};
                    case {'ri','refractiveindex','refractive_index'}
                        nMedium=varargin{ii}{2};
                    otherwise
                        warning('ott:SLM_image_to_bsc:input',['The input is not recognised! The code probably isn''t doing what you want it to as a result!'])
                end
        end
    end
end

NAonm=NA/nMedium;

anglez=asin(NAonm);

%overfit points because I can. This is the angle regridding step.
%everything can be done in one go but it's this way due to legacy code.
nT=min(size(E_ff))*2;
nP=min(size(E_ff))*2;

[theta,phi]=meshgrid(linspace(0,anglez,nT),linspace(0,2*pi,nP));
if NAonm<0
    [theta,phi]=meshgrid(linspace(pi-abs(anglez),pi,nT),linspace(0,2*pi,nP));
end

% tan theta scaling, thin lens appropriate:
wscaling=1/tan(abs(anglez));

Xt=tan(theta).*cos(phi);
Yt=tan(theta).*sin(phi);

if function_theta==1
    % theta scaling:
    wscaling=1/(abs(anglez));
    
    Xt=(theta).*cos(phi);
    Yt=(theta).*sin(phi);
    
end
if function_theta==2
    %sin theta scaling:
    wscaling=1/sin(abs(anglez));
    
    Xt=sin(theta).*cos(phi);
    Yt=sin(theta).*sin(phi);
    
end

%Cartesean coordinates for the paraxial plane. Normalise to BFP:
mXY=min(size(E_ff));
[X,Y]=meshgrid(linspace(-1,1,size(E_ff,2))*size(E_ff,2)/mXY/wscaling*(1+1e-12),linspace(-1,1,size(E_ff,1))*size(E_ff,1)/mXY/wscaling*(1+1e-12));

Exy=interp2(X,Y,E_ff,Xt,Yt);
Exy(isnan(Exy))=0;

% for pointmatching we need the full 4*pi steradians, rescale again:
if ~exist('ntheta','var')
    [Theta,Phi]=angulargrid(2*(Nmax+1),2*(Nmax+1));
else
    [Theta,Phi]=angulargrid(ntheta,nphi);
end

Exy_toolbox=interp2(theta*(1+1e-8),phi,Exy,reshape(Theta,[2*(Nmax+1),2*(Nmax+1)]),reshape(Phi,[2*(Nmax+1),2*(Nmax+1)]));

theta=Theta;
phi=Phi;

Exy_toolbox(isnan(Exy_toolbox))=0;

Ex_toolbox=polarisation(1)*Exy_toolbox(:);
Ey_toolbox=polarisation(2)*Exy_toolbox(:);

if verbose
    figure(1)
    imagesc(abs(Exy_toolbox))
    figure(2)
    plot(theta(:),abs(Ex_toolbox(:)).'/max(abs(Ex_toolbox(:)).'),'o-');
    axis square
    pause
end

if any([ra,az])
    Et=sign(cos(theta)).*Ex_toolbox;
    Ep=Ey_toolbox;
else
    Et=(sign(cos(theta)).*cos(phi).*Ex_toolbox+sign(cos(theta)).*sin(phi).*Ey_toolbox);
    Ep=(-sin(phi).*Ex_toolbox+cos(phi).*Ey_toolbox);
end

e_field=[Et(:);Ep(:)];

mode_indexes=[1:Nmax*(Nmax+2)].';

[nn,mm]=combined_index(mode_indexes);

coefficient_matrix = zeros(length(e_field),2*length(nn));

for n = 1:max(nn)
    ci=find(nn==n);
    
    [~,dtY,dpY]= spharm(n,mm(ci),theta,phi);
    
    coefficient_matrix(:,ci) = [dpY;-dtY] * 1i^(n+1)/sqrt(n*(n+1));
    coefficient_matrix(:,ci+length(nn)) = [dtY;dpY] * 1i^(n)/sqrt(n*(n+1));
    
end

expansion_coefficients = coefficient_matrix \ e_field;

a = expansion_coefficients(1:end/2,:);
b = expansion_coefficients(1+end/2:end,:);

if nargout==2
    nn=a;
    mm=b;
end
if nargout==4
    ci=combined_index(nn,mm);
    a=a(ci);
    b=b(ci);
end
