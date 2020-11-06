classdef PmParaxial < ott.beam.BscInfinite ...
    & ott.beam.properties.Mapping ...
    & ott.beam.properties.Polarisation ...
    & ott.beam.properties.Profile
% Construct a beam using paraxial far-field point matching.
% Inherits from :class:`BscInfinite` and :class:`+properties.Mapping`,
% :class:`+properties.Polarisation` and :class:`+properties.Profile`.
%
% Properties
%   - mapping     -- Paraxial to far-field mapping
%   - polbasis    -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield    -- (2 numeric) Field in theta/phi or x/y directions
%   - basis       -- PM basis ('lg', 'bessel', 'plane', 'auto', or 'full')
%   - truncation_angle -- Maximum angle for beam.  Default: []
%   - data        -- Internal BSC data
%
% Methods
%   - getData         -- Get data for specific Nmax
%   - recalculate     -- Update the internal data for new Nmax
%
% Static methods
%   - InterpProfile   -- Generate a beam profile using interpolation

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    basis
    truncation_angle
  end

  properties (Hidden, SetAccess=protected)
    basisInternal
    truncation_angleInternal
  end

  methods (Static)
    function beam = InterpProfile(X, Y, A, varargin)
      % Construct a beam using interpolation
      %
      % Usage
      %   beam = InterpProfile(X, Y, A, ...)
      %
      % Parameters
      %   - X, Y (NxM numeric) -- Coordinate arrays.  Should be in
      %     ``ndgrid`` format (uses ``griddedInterpolant``).
      %
      %   - A (NxM numeric) -- Amplitude array.
      %
      % All other parameters passed to constructor.

      F = griddedInterpolant(X, Y, A);
      beam = ott.beam.PmParaxial(@(x, y) F(x, y), varargin{:});
    end
  end

  methods
    function beam = PmParaxial(varargin)
      % Construct beam using paraxial point matching
      %
      % Usage
      %   beam = PmParaxial(profile, ...)
      %
      % Optional named parameters
      %   - profile ('uniform' | function_handle | beam) -- Profile for
      %     paraxial point matching.  Default: ``uniform``.
      %
      %   - truncation_angle (numeric) -- Maximum angle for beam.
      %
      %   - basis (enum) -- Basis for pointmatching.  Can be one of
      %     ''lg', 'bessel', 'plane', 'auto', or 'full'.  Default: ``'auto'``.
      %     When 'auto', attempts to choose a basis based on the profile.
      %
      %   - polfield (2 numeric) -- Field in the [theta, phi] or [x, y]
      %     directions.  Default: ``[1, 1i]``.
      %
      %   - polbasis (enum) -- Polarisation basis.  Can be either
      %     'polar' or 'cartesian'.  Default: ``'cartesian'``.
      %
      %   - Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('profile', 'uniform');
      p.addParameter('truncation_angle', []);
      p.addParameter('polfield', [1, 1i]);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('Nmax', 0);
      p.addParameter('basis', 'auto');
      p.parse(varargin{:});

      beam.profile = p.Results.profile;
      beam.truncation_angle = p.Results.truncation_angle;
      beam.polfield = p.Results.polfield;
      beam.polbasis = p.Results.polbasis;
      beam.basis = p.Results.basis;
      beam = beam.recalculate(p.Results.Nmax);
    end
  end

  methods (Hidden)
    function [data, vswfData] = recalculateInternal(beam, Nmax, vswfData)
      % Re-calculate BSC data for specified Nmax

      % TODO
    end
  end

  methods % Getters/setters
    function beam = set.basis(beam, val)
      assert(sum(strcmpi(val, {'lg', 'bessel', 'plane', ...
          'auto', 'full'})) == 1, ...
          'basis must be ''lg'' ''bessel'' ''plane'' ''auto'' or ''full''');
      beam.basisInternal = val;
      beam.data = [];
    end
    function val = get.basis(beam)
      val = beam.basisInternal;
    end

    function beam = set.truncation_angle(beam, val)
      assert(isempty(val) || (isnumeric(val) && isscalar(val)), ...
          'truncation_angle must be numeric scalar');
      beam.truncation_angleInternal = val;
      beam.data = [];
    end
    function val = get.truncation_angle(beam)
      val = beam.truncation_angleInternal;
    end
  end
end

