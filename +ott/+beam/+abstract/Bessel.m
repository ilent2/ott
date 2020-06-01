classdef Bessel < ott.beam.properties.Bessel ...
    & ott.beam.abstract.CastBoth ...
    & ott.beam.properties.InfinitePower
% Abstract description of a Bessel beam.
% Inherits from :class:`ott.beam.properties.Bessel` and :class:`CastBoth`.
%
% Casts
%   Ray
%   vswf.Bsc
%   vswf.Bessel
%   vswf.BesselBasis   -- (Sealed) Cast to array of Bsc

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    angle  % Far-field angle of Bessel beam (radians)
    field  % Field in theta and phi directions
    lmode  % Azimuthal angular momentum number
  end

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

    function beam = ott.beam.Ray(beam, varargin)
      % Construct a Ray beam

      assert(isa(beam, 'ott.beam.abstract.Bessel'), ...
          'First argument must be abstract.Bessel');

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

  methods (Sealed)
    function beam = ott.beam.vswf.BesselBasis(beam, varargin)
      % Cast to vswf.Bessel beam

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be a abstract.Beam');
      ott.utils.nargoutCheck(beam, nargout);

      beam = ott.beam.vswf.BesselBasis.like(beam, varargin{:});
    end
  end

  methods % Getters/setters
    function beam = set.angle(beam, val)
      assert(isnumeric(val) && isscalar(val), 'angle must be numeric scalar');
      beam.angle = val;
    end

    function beam = set.lmode(beam, val)
      assert(isnumeric(val) && isscalar(val), 'lmode must be numeric scalar');
      beam.lmode = val;
    end

    function beam = set.field(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'field must be 2 element numeric');
      beam.field = val(:);
    end
  end
end
