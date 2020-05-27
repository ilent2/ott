classdef PlaneWave < ott.beam.properties.PlaneWaveScalar ...
    & ott.beam.abstract.CastBoth ...
    & ott.beam.properties.InfinitePower
% Abstract representation of a plane wave beam.
% Inherits from :class:`Beam`
% and :class:`ott.beam.properties.PlaneWaveScalar`.
%
% This class creates a single abstract Plane Wave beam.  For more
% efficient arrays of plane wave beams, use :class;`ott.beam.PlaneWave`
% directly.
%
% Supported casts
%   - Beam            -- Default Beam cast, uses PlaneWave
%   - vswf.Bsc        -- Default Bsc cast, uses vswf.PlaneWave
%   - PlaneWave
%   - Ray
%   - vswf.PlaneWave
%   - vswf.PlaneBasis

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = PlaneWave(varargin)
      % Construct a new abstract plane wave representation
      %
      % Usage
      %   beam = PlaneWave(...)
      %
      % Optional named arguments
      %   - position (3 numeric) -- Position of the beam.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Rotation of the beam.
      %     Default: ``eye(3)``.
      %
      %   - field (2 numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %     Default: ``[1, 1i]``.

      beam = beam@ott.beam.properties.PlaneWaveScalar(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to PlaneWave
      beam = ott.beam.PlaneWave(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to vswf.PlaneWave
      beam = ott.beam.vswf.PlaneWave(varargin{:});
    end

    function beam = ott.beam.PlaneWave(beam, varargin)
      % Cast to PlaneWave
      beam = castHelper(@ott.beam.PlaneWave.like, beam, varargin{:});
    end

    function beam = ott.beam.vswf.PlaneWave(beam, varargin)
      % Cast to PlaneWave
      beam = castHelper(@ott.beam.vswf.PlaneWave.like, beam, varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Cast to Ray
      beam = castHelper(@ott.beam.Ray.like, beam, varargin{:});
    end

    function beam = ott.beam.vswf.PlaneBasis(beam, varargin)
      % Cast to PlaneWave
      %
      % Usage
      %   beam = ott.beam.vswf.PlaneBasis(beam, ...)
      %   Additional parameters are passed to vswf.PlaneBasis.like

      beam = castHelper(@ott.beam.vswf.PlaneBasis.like, beam, varargin{:});
    end
  end

  methods (Access=protected)
    function beam = castHelper(cast, beam, varargin)
      % Helper for casts

      assert(isa(beam, 'ott.beam.abstract.PlaneWave'), ...
          'First argument must be a abstract.PlaneWave');
      ott.utils.nargoutCheck(beam, nargout);

      % All other types are PlaneWaveArray's, so add properties
      args = ott.utils.addDefaultParameter('position', [0;0;0], varargin);
      args = ott.utils.addDefaultParameter('rotation', eye(3), args);
      args = ott.utils.addDefaultParameter('origin', [beam.position], args);
      args = ott.utils.addDefaultParameter('field', [beam.field], args);
      args = ott.utils.addDefaultParameter(...
          'directionSet', [beam.rotation], args);

      beam = cast(beam(1), args{:});
    end
  end
end

