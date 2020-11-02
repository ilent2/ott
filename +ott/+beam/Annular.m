classdef Annular < ott.beam.BscInfinite ...
    & ott.beam.properties.Polarisation & ott.beam.properties.Lmode ...
    & ott.beam.properties.Profile
% Construct a VSWF representation of a finite Annular beam
% Inherits from :class:`+ott.+beam.BscInfinite` and
% :class:`+ott.+beam.+properties.Polarisation` and
% :class:`+ott.+beam.+properties.Profile` and
% :class:`+ott.+beam.+properties.Lmode`.
%
% Represents beams internally as an array of Bessel-like beams, hence
% the class inherits from :class:`+ott.+beam.BscInfinite`.
%
% Properties
%   - polbasis      -- Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- Polarisation field (theta/phi or x/y)
%   - theta         -- Low and upper angles of the annular
%   - profile       -- Anular profile (function_handle | beam | 'uniform')
%   - lmode         -- Orbital angular momentum number
%
% Static methods
%   - InterpProfile   -- Generate a beam profile using interpolation

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

      F = griddedInterpolant(theta, A):
      beam = ott.beam.Annular([theta(1), theta(end)], ...
          @(t) F(t), varargin{:});
    end
  end

  methods
    function beam = Annular(varargin)
      % Construct a new VSWF representation of a Annular beam.
      %
      % Usage
      %   beam = Annular(theta, profile, ...)
      %
      % Optional named parameters
      %   - theta (2 numeric) -- Angles specifying range for annular (radians).
      %     Default: ``[2*pi/8, 3*pi/8]``
      %
      %   - profile ('uniform' | function_handle | beam) -- Specifies the
      %     angular profile of the beam.  Can be either a
      %     :class:`Beam` or :class:`Bsc` instance, a function handle
      %     with the signature ``@(theta) intensity`` or ``'uniform'``.
      %     Default: ``'uniform'``.
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
      %   - initial_Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('theta', [2*pi/8, 3*pi/8]);
      p.addOptional('profile', 'uniform');
      p.addParameter('lmode', 0);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('polfield', [1, 1i]);
      p.addParameter('initial_Nmax', 0);
      p.parse(varargin{:});

      beam.theta = p.Results.theta;
      beam.profile = p.Results.profile;
      beam.lmode = p.Results.lmode;
      beam.polbasis = p.Results.polbasis;
      beam.polfield = p.Results.polfield;

      beam = beam.recalculate(p.Results.initial_Nmax);
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
