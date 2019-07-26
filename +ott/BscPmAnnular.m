classdef BscPmAnnular < ott.BscPointmatch
%BscPmAnnular gennerate a beam with an annular far-field profile
%
% Properties:
%   a               (Bsc) Beam shape coefficients a vector
%   b               (Bsc) Beam shape coefficients b vector
%   type            (Bsc) Beam type (incoming, outgoing or scattered)
%   NA              Numerical aperture of inner and outer edge [r1, r2]
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
% See also BscPmAnnular, ott.Bsc and ott.BscPmGauss
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    NA            % Numerical aperture of inner and outer edge [r1, r2]
  end
  
  methods
    function beam = BscPmAnnular(NA, varargin);
      % BscPmAnnular gennerate BSC with an annular far-field profile
      %
      % beam = BscPmParaxial(NA, ...) generates beam with an annular
      % far-field profile with minimum and maximum annular edge angles
      % specified as a numerical aperture range NA [min, max].
      % Default is a uniform intensity, flat phase profile.
      %
      % Optional named parameters:
      %     profile   opt    Specify the profile of the beam.  This can
      %        be another beam or a 1-D array representing the radial
      %        profile in equal angle steps from NAmin to NAmax.
      %     grid  [ntheta, nphi] Number of grid points
      %        to match far-field and BSC values at (default: [])
      %     oam       num    Orbital angular momentum number (default: 0)
      %        Only used with 1-D and oam beam profiles.
      %
      %     verbose   bool   Display additional output (default: false)
      %     Nmax      num    Truncation number for BSC (default: 30)
      %     polarisation  [x,y]  Polarisation of far-field (default: [1, 0])
      %     omega     num    Angular frequency of beam (default: 2*pi)
      %     wavelength0 num  Wavelength of beam in vacuum (default: 1)
      %     k_medium    num  Wave-number of beam in medium (default: [])
      %     index_medium num Refractive index of medium (default: [])
      %     wavelength_medium num Wavelength of beam in medium (default: [])
      %
      %     keep_coefficient_matrix bool True to calculate and keep the
      %       inverted coefficient matrix. (default: false)
      %     invert_coefficient_matrix bool True to use the inverted
      %       coefficient matrix even if mldivide is available.
      %       (default: keep_coefficient_matrix).

      % Based on BscPmParaxial
      
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
      p.addParameter('profile', []);
      p.addParameter('oam', 0);
      p.addParameter('grid', []);
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

      % Update progress callback
      p.Results.progress_callback(0.1);

      NAonm=-NA/nMedium;

      anglez=asin(NAonm);
      
      % In order to deal with 0 we need same hemisphere
      % Otherwise we need to change to theta instead of NA as input
      assert(prod(sign(anglez)) >= 0, ...
          'NA must be both in same hemisphere');
      if any(anglez < 0)
        anglez(anglez == 0) = pi;
      end
      
      % Negative NAonm refers to oposite hemisphere
      anglez(NAonm < 0) = pi - abs(anglez(NAonm < 0));
      assert(numel(anglez) == 2, 'NA must be 2 values');
      assert(min(anglez) >= 0 && max(anglez) <= pi, ...
          'anglez outside range [0,pi]');
      
      % Calculate grid for point matching (need 4*pi steradians)
      if isempty(p.Results.grid)
          [theta,phi]=ott.utils.angulargrid(2*(Nmax+1),2*(Nmax+1));
      else
          [theta,phi]=ott.utils.angulargrid(p.Results.grid(1),p.Results.grid(2));
      end
      
      % Generate profile (e_field)
      if isempty(p.Results.profile)
        
        % Profile is uniform (posibly with oam)
        profile = exp(2*pi*1i*p.Results.oam*phi);
        profile(theta < min(anglez)) = 0;
        profile(theta > max(anglez)) = 0;
        
        % Convert scalar to x/y polarised far-field
        Ex_toolbox=polarisation(1)*profile(:);
        Ey_toolbox=polarisation(2)*profile(:);
        
        % Convert to spherical coordinates
        Et=(sign(cos(theta)).*cos(phi).*Ex_toolbox+sign(cos(theta)).*sin(phi).*Ey_toolbox);
        Ep=(-sin(phi).*Ex_toolbox+cos(phi).*Ey_toolbox);

        % Generate e_field for point matching
        e_field=[Et(:);Ep(:)];
        
      elseif isa(p.Results.profile, 'ott.Bsc')
        
        % Beam object supplied (no oam)
        [e_field, ~] = p.Results.profile.farfield(theta, phi);
        
        % Remove parts outside annular
        e_field(:, theta < min(anglez)) = 0;
        e_field(:, theta > max(anglez)) = 0;
        
        % Strip R component and package
        e_field = [e_field(2, :).'; e_field(3, :).'];
        
      elseif min(size(p.Results.profile)) == 1
        
        % Radial profile of beam (posibly with oam)
        
        % Start with a uniform (with oam) profile
        profile = exp(2*pi*1i*p.Results.oam*phi);
        profile(theta < min(anglez)) = 0;
        profile(theta > max(anglez)) = 0;
        
        % Add the modulation
        profileAngles = linspace(min(anglez), max(anglez), numel(p.Results.profile));
        profileMod = interp1(profileAngles, p.Results.profile, theta);
        profileMod(isnan(profileMod)) = 0;
        profile = profile .* profileMod;
        
        % Convert scalar to x/y polarised far-field
        Ex_toolbox=polarisation(1)*profile(:);
        Ey_toolbox=polarisation(2)*profile(:);
        
        % Convert to spherical coordinates
        Et=(sign(cos(theta)).*cos(phi).*Ex_toolbox+sign(cos(theta)).*sin(phi).*Ey_toolbox);
        Ep=(-sin(phi).*Ex_toolbox+cos(phi).*Ey_toolbox);

        % Generate e_field for point matching
        e_field=[Et(:);Ep(:)];
        
      else
        
        error('Unknown profile argument value');
        
      end

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

