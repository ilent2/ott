classdef Annular < ott.beam.BscWBessel ...
    & ott.beam.properties.Polarisation & ott.beam.properties.Lmode ...
    & ott.beam.properties.Profile
% Construct a VSWF representation of a finite Annular beam
% Inherits from :class:`+ott.+beam.BscWBessel` and
% :class:`+ott.+beam.+properties.Polarisation` and
% :class:`+ott.+beam.+properties.Profile` and
% :class:`+ott.+beam.+properties.Lmode`.
%
% Represents beams internally as an array of Bessel-like beams.
% Applies the beam weights when the data is requested.
%
% Properties
%   - polbasis      -- Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- Polarisation field (theta/phi or x/y)
%   - theta         -- Low and upper angles of the annular
%   - profile       -- Anular profile (function_handle)
%   - lmode         -- Orbital angular momentum number
%
% Static methods
%   - InterpProfile   -- Generate a beam profile using interpolation
%   - BeamProfile     -- Generate a beam profile from another beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    theta         % Low and upper angles of the annular
  end

  properties (Hidden, SetAccess=protected)
    thetaInternal
  end

  methods (Static)
    function beam = InterpProfile(theta, A, varargin)
      % Construct a beam using interpolation
      %
      % Usage
      %   beam = InterpProfile(theta, A, ...)
      %
      % Parameters
      %   - theta (N numeric) -- Angular coordinates (radians).
      %
      %   - A (N numeric) -- Amplitude array.
      %
      % All other parameters passed to constructor.

      F = griddedInterpolant(theta, A);
      profile = @(t) F(t);
      beam = ott.beam.Annular([theta(1), theta(end)], ...
          'profile', profile, varargin{:});
    end

    function beam = BeamProfile(theta, beam, varargin)
      % Construct a annular beam masking another beam.
      %
      % Usage
      %   beam = BeamProfile(theta, other_beam, ...)
      %
      % Parameters
      %   - theta (2 numeric) -- Theta range for annular.
      %
      %   - other_beam (ott.beam.Beam) -- The beam to use.  Uses the
      %     beams efarfield method to calculate fields.
      %
      % All other parameters passed to constructor.

      profile = @(theta) vecnorm(beam.efarfield(...
          [theta(:), 0*theta(:)].').vxyz);

      beam = ott.beam.Annular(theta, 'profile', profile, varargin{:});
    end
  end

  methods
    function beam = Annular(varargin)
      % Construct a new VSWF representation of a Annular beam.
      %
      % Usage
      %   beam = Annular(theta, ...)
      %
      % Optional named parameters
      %   - theta (2 numeric) -- Angles specifying range for annular (radians).
      %     Default: ``[2*pi/8, 3*pi/8]``
      %
      %   - profile (function_handle) -- Specifies the angular profile
      %     of the beam.  Default: ``@(theta) ones(size(theta))``.
      %
      %   - lmode (numeric) -- Azimuthal mode number.  Adds an orbital
      %     component to the beam.  Default: ``0`` (no OAM).
      %
      %   - polbasis (enum) -- Polarisation basis.  Either 'polar' or
      %     'cartesian'.  Default: ``'cartesian'``.
      %
      %   - polfield (2 numeric) -- Field in either the theta/phi or
      %     x/y directions.  Default: ``[1, 1i]``.
      %
      %   - Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('theta', [2*pi/8, 3*pi/8], @isnumeric);
      p.addParameter('profile', @(val) ones(size(val)));
      p.addParameter('lmode', 0);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('polfield', [1, 1i]);
      p.addParameter('Nmax', 0);
      p.parse(varargin{:});

      beam.theta = p.Results.theta;
      beam.profile = p.Results.profile;
      beam.lmode = p.Results.lmode;
      beam.polbasis = p.Results.polbasis;
      beam.polfield = p.Results.polfield;

      beam = beam.recalculate(p.Results.Nmax);
    end

    function beam = recalculate(beam, Nmax)
      % Re-calculate BSC data for specified Nmax

      % Get angles for Bessel beams
      nTheta = 2*(Nmax+1);
      otheta = linspace(beam.theta(1), beam.theta(2), nTheta+1);
      otheta = otheta(1:end-1) + diff(otheta(1:2))./2;

      % Get weights of each segment
      beam.besselWeights = beam.profile(otheta) ...
          .* beam.wavenumber .* sin(otheta) ./ nTheta;

      % TODO: Different polarisation support

      % Calculate Bsc data
      beam.data = ott.bsc.Annular.FromBessel(Nmax, otheta, ...
          beam.polfield(:), beam.lmode);
    end
  end

  methods % Getters/setters
    function beam = set.theta(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'theta must be 2 element numeric');
      beam.thetaInternal = val;
      beam.data = [];
    end
    function val = get.theta(beam)
      val = beam.thetaInternal;
    end
  end
end
