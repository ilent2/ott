classdef BscPmParaxial < ott.BscPointmatch
%BscPmParaxial calculate representation from farfield/paraxial beam
%
% BscPmParaxial properties:
%
% BscPmParaxial methods:
%
% Based on paraxial_to_bsc from ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
  end

  methods
    function beam = BscPmParaxial(NA, E_ff, varargin);
      % PARAXIAL_TO_BSC Takes complex amplitudes in a plane (rho \propto \theta)
      % and maps it onto an angular grid with extent dermined by NA, and solves 
      % the BSCs.
      % NOTE: Default behaviour is to have the BFP radius \propto sin(theta).
      %
      % TODO: Documentation
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
      % NOTE: This current version will best work for "clean" beam modes, that
      % is, the desired field AFTER spatial filtering (If modelling an SLM/DMD).
      % In addition, a standard G&L algorithm will produce abberations in the 
      % output BSC.

      p = inputParser;
      p.addParameter('verbose', false);
      p.addParameter('Nmax', 30);
      p.addParameter('index_medium', 1.33);
      p.addParameter('polarisation', [ 1 0 ]);
      p.addParameter('grid', []);
      p.addParameter('mapping', 'sintheta');
      p.parse(varargin{:});

      verbose = p.Results.verbose;
      Nmax = p.Results.Nmax;
      nMedium = p.Results.index_medium;

      polarisation = p.Results.polarisation;

      % TODO: These parameters were set from input
      ra=false;
      az=false;

      function_theta=[];
      switch lower(p.Results.mapping)
        case 'sintheta'
          function_theta=2;
        case 'tantheta'
          function_theta=0;
        case 'theta'
          function_theta=1;
        otherwise
          error('ott:BscPmParaxial:mapping', 'Unrecognized mapping value');
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
      if isempty(p.Results.grid)
          [Theta,Phi]=ott.utils.angulargrid(2*(Nmax+1),2*(Nmax+1));
      else
          [Theta,Phi]=ott.utils.angulargrid(grid(1),grid(2));
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

      [nn,mm]=ott.utils.combined_index(mode_indexes);

      % Do the point matching and store the result
      [beam.a, beam.b] = beam.bsc_farfield(nn, mm, e_field, theta, phi);

      beam.type = 'incomming';
      beam.k_medium = 2*pi*nMedium;
    end
  end
end

