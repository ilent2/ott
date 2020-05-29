classdef Bessel < ott.beam.properties.Bessel ...
    & ott.beam.abstract.CastBoth ...
    & ott.beam.properties.InfinitePower
% Abstract description of a Bessel beam.
% Inherits from :class:`ott.beam.properties.Bessel` and :class:`CastBoth`.
%
% Casts
%   Ray
%   vswf.Bsc           -- Default Bsc cast, uses vswf.Bessel
%   vswf.Bessel
%   vswf.BesselBasis

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Bessel(varargin)
      % Construct a new abstract Bessel beam
      %
      % Usage
      %   beam = Bessel(theta, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - theta (numeric) -- Angle of the Bessel beam from beam axis.

      beam = beam@ott.beam.properties.Bessel(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to vswf.Bessel beam
      beam = ott.beam.vswf.Bessel(varargin{:});
    end

    function beam = ott.beam.vswf.Bessel(varargin)
      % Cast to vswf.Bessel beam
      beam = castHelper(@ott.beam.vswf.Bessel.like, varargin{:});
    end

    function beam = ott.beam.vswf.BesselBasis(varargin)
      % Cast to vswf.Bessel beam
      beam = castArrayHelper(@ott.beam.vswf.BesselBasis.like, varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Construct a Ray beam

      assert(isa(beam, 'ott.beam.abstract.Bessel'), ...
          'First argument must be abstract.Bessel');

      if numel(beam) > 1
        oldbeam = beam;
        beam = ott.beam.Array('coherent', size(beam));
        for ii = 1:numel(oldbeam)
          beam(ii) = ott.beam.Ray(oldbeam(ii), varargin{:});
        end
      else
        phi = linspace(0, 2*pi, 100);
        theta = beam.angle.*ones(size(phi));
        radius = 10.*ones(size(phi));

        rtp = [radius(:), theta(:), phi(:)].';
        vrtp = [0*radius(:), ones(size(theta(:))), 0*phi(:)].';
        [vxyz, xyz] = ott.utils.rtpv2xyzv(vrtp, rtp);

        Etp = [1;0].*ones(1, size(xyz, 2));

        directionSet = ott.beam.Ray.DirectionSet(-xyz, vxyz);
        beam = ott.beam.Ray.like(beam, 'origin', xyz, ...
            'directionSet', directionSet, 'field', Etp, varargin{:});
      end
    end
  end
end
