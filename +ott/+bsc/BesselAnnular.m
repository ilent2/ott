classdef BscBesselAnnular < ott.Bsc
%BscBesselAnnular gennerate a beam with an annular far-field profile
% This code uses annular beams and may be significantly faster than PM.
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
% See also BscBesselAnnular, and ott.BscPmAnnular
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    NA            % Numerical aperture of inner and outer edge [r1, r2]
  end
  
  methods
    function beam = BscBesselAnnular(NA, varargin);
      % BscPmAnnular gennerate BSC with an annular far-field profile
      %
      % beam = BscPmParaxial(NA, ...) generates beam with an annular
      % far-field profile with minimum and maximum annular edge angles
      % specified as a numerical aperture range NA [min, max].
      % Default is a uniform intensity, flat phase profile.
      %
      % Optional named parameters:
      %     profile   opt    Specify the profile of the beam.  This is
      %        1-D array representing the radial profile in equal
      %        angle steps from NAmin to NAmax.
      %     ntheta    num    Number of angular points to match
      %        to match far-field and BSC values at (default: 100)
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

      % Based on BscPmParaxial
      
      p = inputParser;
      p.addParameter('verbose', false);
      p.addParameter('progress_callback', @(x) []);
      p.addParameter('Nmax', 30);

      p.addParameter('omega', 2*pi);
      p.addParameter('wavelength0', 1);
      p.addParameter('k_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength_medium', []);

      p.addParameter('polarisation', [ 1 0 ]);
      p.addParameter('profile', []);
      p.addParameter('oam', 0);
      p.addParameter('ntheta', 100);
      p.parse(varargin{:});

      verbose = p.Results.verbose;
      Nmax = p.Results.Nmax;
      beam.k_medium = ott.Bsc.parser_k_medium(p, 2*pi);
      beam.omega = p.Results.omega;

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
      
      % Calculate angles for bessel beams
      theta = linspace(min(anglez), max(anglez), p.Results.ntheta);
      
      % Generate bessel beams
      beams = ott.BscBessel(Nmax, theta, ...
          'k_medium', beam.k_medium, ...
          'polarisation', polarisation, ...
          'lmode', p.Results.oam);
      
      % Generate profile (e_field)
      if isempty(p.Results.profile)
        
        area = sin(theta);
        beam.a = beams.a * area(:);
        beam.b = beams.b * area(:);
        
      elseif min(size(p.Results.profile)) == 1
        
        % Radial profile of beam (posibly with oam)
        
        error('Not yet implemented');
        
      else
        
        error('Unknown profile argument value');
        
      end

      beam.type = 'incident';
      beam.basis = 'regular';
      
      % Update progress callback
      p.Results.progress_callback(1.0);
    end
  end
end

