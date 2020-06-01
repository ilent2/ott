classdef TopHat < ott.beam.properties.TopHat ...
    & ott.beam.abstract.CastBoth
% Abstract description of Top-hat beams
% Inherits from :class:`ott.beam.properties.Bessel` and :class:`CastBoth`.
%
% This class describes a collimated Top-Hat beam.
%
% Casts
%   - Beam      -- Cast to abstract.MaskedNearfield3d
%   - vswf.Pointmatch -- Uses vswf.NearfieldPm
%   - vswf.NearfieldPm -- Uses abstract.MaskedNearfield3d
%   - Ray
%   - abstract.MaskedNearfield3d

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    power
  end

  methods
    function beam = TopHat(varargin)
      % Construct a new abstract TopHat beam
      %
      % Usage
      %   beam = TopHat(radius, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - radius (numeric) -- Angle of the Bessel beam from beam axis.
      %
      % Optional named arguments
      %   - polarisation (2 numeric) -- Polarisation of beam.
      %     Default: ``[1; 0]``.

      p = inputParser;
      p.addOptional('radius', [], @isnumeric);
      p.addParameter('polarisation', [1;0]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.TopHat(...
          'radius', p.Results.radius, ...
          'polarisation', p.Results.polarisation, ...
          unmatched{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to MaskedNearfield3d
      beam = ott.beam.abstract.MaskedNearfield3d(varargin{:});
    end

    function beam = ott.beam.vswf.Pointmatch(varargin)
      % Use NearfieldPm
      beam = ott.beam.vswf.NearfieldPm(varargin{:});
    end

    function beam = ott.beam.vswf.NearfieldPm(beam, varargin)
      beam = ott.beam.abstract.MaskedNearfield3d(beam);
      beam = ott.beam.vswf.NearfieldPm(beam, varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Construct a Ray beam

      assert(isa(beam, 'ott.beam.abstract.TopHat'), ...
          'First argument must be abstract.TopHat');

      radius = linspace(0, beam.radius, 100);
      phi = linspace(0, 2*pi, 100);
      [R, P, Z] = meshgrid(radius, phi, 0);
      X = R .* sin(P);
      Y = R .* cos(P);
      xyz = [X(:), Y(:), Z(:)].';

      vdir = [0;0;1].*ones(1, size(xyz, 2));
      vpol = [1;0;0].*ones(1, size(xyz, 2));

      directionSet = ott.beam.Ray.DirectionSet(vdir, vpol);
      beam = ott.beam.Ray.like(beam, 'origin', xyz, ...
          'directionSet', directionSet, ...
          'field', beam.polarisation, varargin{:});
    end

    function beam = ott.beam.abstract.MaskedNearfield3d(beam, varargin)
      % Construct a masked near-field instance

      assert(isa(beam, 'ott.beam.abstract.TopHat'), ...
          'First argument must be abstract.TopHat');

      plane = ott.beam.abstract.PlaneWave.FromDirection(...
          'origin', [0;0;0], 'direction', [0;0;1], ...
          'polarisation', [1;0;0], 'field', beam.polarisation);

      mask = @(xyz) vecnorm(xyz(1:2, :)) <= beam.radius;

      beam = ott.beam.abstract.MaskedNearfield3d(mask, plane);
    end
  end

  methods % Getters/setters
    function power = get.power(beam)
      power = pi.*beam.radius.^2 .* vecnorm(beam.polarisation);
    end
    function beam = set.power(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'power must be numeric scalar');
      beam.polarisation = beam.polarisation ...
          ./ vecnorm(bean.polarisation) .* val ./ (pi*beam.radius.^2);
    end
  end
end

