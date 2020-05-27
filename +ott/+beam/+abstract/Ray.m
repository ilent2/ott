classdef Ray < ott.beam.properties.PlaneWaveScalar ...
    & ott.beam.abstract.Beam
% Specialisation describing abstract geometric ray beams.
% Inherits from :class:`ott.beam.properties.Ray` and :class:`Beam`.
%
% This class is only provided for consistency.  To create a Ray beam,
% it is recommended to create either a :class:`ott.beam.Ray` directly
% or create another beam and cast to a `ott.beam.Ray`, for example
%
% .. code:: matlab
%   beam = ott.beam.abstract.Gaussian(1.0);
%   rays = ott.beam.Ray(beam);
%
% Supported casts
%   - Beam                -- Default Beam cast, uses Ray
%   - Ray
%   - abstract.PlaneWave

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    power
  end

  methods
    function beam = Ray(varargin)
      % Construct a new abstract single ray representation.
      %
      % Usage
      %   beam = Ray(...)
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
      % Cast to Ray
      beam = ott.beam.Ray(varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Cast to Ray
      beam = castHelperArray(@ott.beam.Ray.like, beam, varargin{:});
    end

    function beam = ott.beam.abstract.PlaneWave(beam, varargin)
      % Cast to abstract.PlaneWave
      beam = castHelper(@ott.beam.abstract.PlaneWave.like, beam, varargin{:});
    end
  end

  methods (Access=protected)
    function beam = castHelperArray(cast, beam, varargin)
      % Helper for casts

      assert(isa(beam, 'ott.beam.abstract.Ray'), ...
          'First argument must be a abstract.Ray');
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

  methods % Getters/setters
    function p = get.power(beam)
      p = beam.intensity;
    end
  end
end
