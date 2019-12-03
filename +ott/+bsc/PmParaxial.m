classdef BscPmParaxial < ott.BscPointmatch
%BscPmParaxial calculate representation from farfield/paraxial beam
%
% Properties:
%   a               (Bsc) Beam shape coefficients a vector
%   b               (Bsc) Beam shape coefficients b vector
%   type            (Bsc) Beam type (incoming, outgoing or scattered)
%
% Methods:
%   translateZ      (Bsc) Translates the beam along the z axis
%   translateXyz    (Bsc) Translation to xyz using rotations and z translations
%   translateRtp    (Bsc) Translation to rtp using rotations and z translations
%   farfield        (Bsc) Calculate fields in farfield
%   emFieldXyz      (Bsc) Calculate fields at specified locations
%   set.Nmax        (Bsc) Resize the beam shape coefficient vectors
%   get.Nmax        (Bsc) Get the current size of the beam shape vectors
%   getCoefficients (Bsc) Get the beam coefficients [a, b]
%   getModeIndices  (Bsc) Get the mode indices [n, m]
%   power           (Bsc) Calculate the power of the beam
%
% Based on paraxial_to_bsc from ottv1.
%
% See also BscPmParaxial, ott.Bsc and examples/slm_to_focalplane
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
  end

  methods
    function beam = BscPmParaxial(NA, E_ff, varargin);
      % BscPmParaxial generates BSC from far-field complex field amplitudes
      %
      % beam = BscPmParaxial(NA, Eff, ...) generates beam with specified
      % numerical aperture (NA) from the complex E-field in
      % the paraxial/far-field (Eff).
      %
      % Optional named parameters
      %   - verbose   bool   Display additional output (default: false)
      %   - Nmax      num    Truncation number for BSC (default: 30)
      %
      %   - polarisation  [x,y]  Polarisation of far-field (default: [1, 0])
      %     Ignored if Eff is a MxNx2 matrix.
      %
      %   - radius (numeric) -- Back aperture radius (pixels).
      %     Default: ``min([size(E_ff, 1), size(E_ff, 2)])/2``.
      %
      %   - grid  [ntheta, nphi] Number of grid points for
      %     to match far-field and BSC values at (default: [])
      %
      %   - mapping (enum) -- Determines how Eff is mapped to back plane of
      %     objective.  Options are:
      %    -  'sintheta' (default)  radius \propto \sin(\theta)
      %    -  'tantheta'            radius \propto \tan(\theta)
      %    -  'theta'               radius \propto \theta
      %
      %   - omega     num    Angular frequency of beam (default: 2*pi)
      %     wavelength0 num  Wavelength of beam in vacuum (default: 1)
      %     k_medium    num  Wave-number of beam in medium (default: [])
      %     index_medium num Refractive index of medium (default: [])
      %     wavelength_medium num Wavelength of beam in medium (default: [])
      %
      %   - keep_coefficient_matrix bool True to calculate and keep the
      %     inverted coefficient matrix. (default: false)
      %   - invert_coefficient_matrix bool True to use the inverted
      %     coefficient matrix even if mldivide is available.
      %     (default: keep_coefficient_matrix).
      %
      % NOTE: This current version will best work for "clean" beam modes, that
      % is, the desired field AFTER spatial filtering (If modelling an SLM/DMD).
      % In addition, a standard G&L algorithm will produce abberations in the 
      % output BSC.

      p = inputParser;
      p.addParameter('verbose', false);
      p.addParameter('progress_callback', @(x) []);
      p.addParameter('beamData', []);
      p.addParameter('keep_coefficient_matrix', false);
      p.addParameter('invert_coefficient_matrix', []);    % Default arg bellow
      p.addParameter('Nmax', 30);

      p.addParameter('omega', 2*pi);
      p.addParameter('wavelength0', 1);
      p.addParameter('k_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength_medium', []);

      p.addParameter('polarisation', [ 1 0 ]);
      p.addParameter('radius', []);
      p.addParameter('grid', []);
      p.addParameter('mapping', 'sintheta');
      p.parse(varargin{:});

      verbose = p.Results.verbose;
      Nmax = p.Results.Nmax;
      beam.k_medium = ott.Bsc.parser_k_medium(p, 2*pi);
      beam.omega = p.Results.omega;

      % Handle default argument for invert_coefficient_matrix
      invert_coefficient_matrix = p.Results.invert_coefficient_matrix;
      if isempty(invert_coefficient_matrix)
        invert_coefficient_matrix = p.Results.keep_coefficient_matrix;
      end

      if isempty(p.Results.index_medium)
        nMedium = 1.0;
      else
        nMedium = p.Results.index_medium;
      end

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
      
      % Update progress callback
      p.Results.progress_callback(0.1);

      NAonm=NA/nMedium;

      anglez=asin(NAonm);
      
      assert(ndims(E_ff) == 3 || ismatrix(E_ff), 'E_ff must be 2 or 3 dimensional array');
      assert(size(E_ff, 3) == 1 || size(E_ff, 3) == 2, 'E_ff must be NxM or NxMx2 matrix');

      %overfit points because I can. This is the angle regridding step.
      %everything can be done in one go but it's this way due to legacy code.
      nT=min([size(E_ff, 1), size(E_ff, 2)])*2;
      nP=min([size(E_ff, 1), size(E_ff, 2)])*2;

      [theta,phi]=meshgrid(linspace(0,anglez,nT),linspace(0,2*pi,nP));
      if NAonm<0
          [theta,phi]=meshgrid(linspace(pi-abs(anglez),pi,nT),linspace(0,2*pi,nP));
      end

      % tan theta scaling, thin lens appropriate:
      wscaling=1/tan(abs(anglez));

      Xt=tan(theta).*cos(phi);
      Yt=tan(theta).*sin(phi);

      if function_theta==2
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
      if isempty(p.Results.radius)
        mXY=min([size(E_ff, 1), size(E_ff, 2)]);
      else
        mXY = 2*p.Results.radius;
      end
      [X,Y]=meshgrid(linspace(-1,1,size(E_ff,2))*size(E_ff,2)/mXY/wscaling*(1+1e-12),...
        linspace(-1,1,size(E_ff,1))*size(E_ff,1)/mXY/wscaling*(1+1e-12));

      Exy = zeros([size(Xt), size(E_ff, 3)]);
      for ii = 1:size(E_ff, 3)
        Exy(:, :, ii) = interp2(X,Y,E_ff(:, :, ii),Xt,Yt);
      end
      Exy(isnan(Exy))=0;

      % for pointmatching we need the full 4*pi steradians, rescale again:
      if isempty(p.Results.grid)
          [Theta,Phi]=ott.utils.angulargrid(2*(Nmax+1),2*(Nmax+1));
      else
          [Theta,Phi]=ott.utils.angulargrid(grid(1),grid(2));
      end

      Exy_toolbox = zeros([2*(Nmax+1),2*(Nmax+1)]);
      for ii = 1:size(E_ff, 3)
        Exy_toolbox(:, :, ii)=interp2(theta*(1+1e-8),phi,Exy(:, :, ii),...
          reshape(Theta,[2*(Nmax+1),2*(Nmax+1)]),reshape(Phi,[2*(Nmax+1),2*(Nmax+1)]));
      end
      
      % Update progress callback
      p.Results.progress_callback(0.2);

      theta=Theta;
      phi=Phi;

      Exy_toolbox(isnan(Exy_toolbox))=0;

      if size(E_ff, 3) == 2
        Ex_toolbox = reshape(Exy_toolbox(:, :, 1), [], 1);
        Ey_toolbox = reshape(Exy_toolbox(:, :, 2), [], 1);
      else
        Ex_toolbox=polarisation(1)*Exy_toolbox(:);
        Ey_toolbox=polarisation(2)*Exy_toolbox(:);
      end

      if verbose
          figure(1)
          imagesc(abs(Exy_toolbox(:, :, 1)));
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

      % Get a previous coefficient matrix
      if ~isempty(p.Results.beamData)
        if isa(p.Results.beamData, 'ott.BscPmParaxial')

          % Check properties of beams match
          % TODO: Check e_field/theta/phi locations
          obeam = p.Results.beamData;
          assert(Nmax == obeam.Nmax, 'Nmax of beams must match');

          icm = p.Results.beamData.inv_coefficient_matrix;
          if isempty(icm)
            warning('ott:BscPmParaxial:BscPmParaxial:empty_icm', ...
                'Not inverse coefficient matrix data associated with beam');
          end
        else
          icm = p.Results.beamData;
        end
      else
        icm = [];
      end

      % Do the point matching and store the result
      if p.Results.keep_coefficient_matrix
        [beam.a, beam.b, icm] = beam.bsc_farfield(nn, mm, e_field, ...
            theta, phi, ...
            'invert_coefficient_matrix', invert_coefficient_matrix, ...
            'inv_coefficient_matrix', icm);
        beam.inv_coefficient_matrix = icm;
      else
        [beam.a, beam.b] = beam.bsc_farfield(nn, mm, e_field, ...
            theta, phi, ...
            'invert_coefficient_matrix', invert_coefficient_matrix, ...
            'inv_coefficient_matrix', icm);
        beam.inv_coefficient_matrix = [];
      end

      beam.type = 'incident';
      beam.basis = 'regular';
      
      % Update progress callback
      p.Results.progress_callback(1.0);
    end
  end
end

