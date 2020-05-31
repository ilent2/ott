classdef (Abstract) PlaneWaveStem < ott.beam.properties.PlaneWave ...
    & ott.beam.abstract.CastBoth ...
% Abstract base class for representation of a plane wave-like beams.
% Inherits from :class:`CastBoth`
% and :class:`ott.beam.properties.PlaneWave`.
%
% Properties
%   - field           -- Field in polarisation directions
%   - origin          -- Alias for position
%   - directionSet    -- Alias for rotation
%
% Supported casts
%   - Beam            -- Raises an error, should be overloaded in sub-class
%   - PlaneWave       -- (Sealed)
%   - Ray             -- (Sealed)
%   - abstract.Ray
%   - abstract.PlaneWave
%   - vswf.Bsc        -- Uses vswf.PlaneWave
%   - vswf.PlaneWave
%   - vswf.PlaneBasis -- (Sealed)
%   - Negative        -- (Sealed) Negates fields

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent, SetAccess=protected)
    directionSet   % Set of direction vectors describing orientation
    origin         % Position used to calculate beam phase offset
  end

  properties (SetAccess=protected)
    field          % Field in polarisation directions
  end

  methods
    function beam = PlaneWaveStem(varargin)
      % Construct a new abstract plane wave or ray representation
      %
      % Usage
      %   beam = PlaneWaveStem(...)
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
      %     Default: ``[1; 1i]``.
      %
      %   - directionSet -- Alias for rotation.
      %   - origin -- Alias for position.
      %
      % See also :meth:`FromDirection` and :meth:`DirectionSet`.

      p = inputParser;
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.addParameter('origin', []);
      p.addParameter('directionSet', []);
      p.addParameter('field', [1; 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get origin from inputs
      if ~any(strcmpi(p.UsingDefaults, 'origin')) ...
          && ~any(strcmpi(p.UsingDefaults, 'position'))
        warning('ott:beam:abstract:PlaneWave:position_and_origin', ...
            'Both position and origin specified, ignoring position');
        origin = p.Results.origin;
      elseif ~any(strcmpi(p.UsingDefaults, 'origin'))
        origin = p.Results.origin;
      else
        origin = p.Results.position;
      end

      % Get directionSet from inputs
      if ~any(strcmpi(p.UsingDefaults, 'directionSet')) ...
          && ~any(strcmpi(p.UsingDefaults, 'rotation'))
        warning('ott:beam:abstract:PlaneWaveStem:directionSet_rotation', ...
            'Both directionSet and rotation specified, ignoring rotation');
        directionSet = p.Results.directionSet;
      elseif ~any(strcmpi(p.UsingDefaults, 'directionSet'))
        directionSet = p.Results.directionSet;
      else
        directionSet = p.Results.rotation;
      end

      beam = beam@ott.beam.properties.PlaneWave(unmatched{:}, ...
          'directionSet', directionSet, 'origin', origin, ...
          'field', p.Results.field);
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to PlaneWave
      error('Method should be overloaded in sub-class');
    end

    %
    % Abstract
    %

    function beam = ott.beam.abstract.Ray(varargin)
      % Cast to abstract ray
      beam = castHelper(@ott.beam.abstract.Ray.like, varargin{:});
    end

    function beam = ott.beam.abstract.PlaneWave(varargin)
      % Cast to abstract ray
      beam = castHelper(@ott.beam.abstract.PlaneWave.like, ...
          varargin{:});
    end

    %
    % VSWF
    %

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to vswf.PlaneWave
      beam = ott.beam.vswf.PlaneWave(varargin{:});
    end

    function beam = ott.beam.vswf.PlaneWave(beam, varargin)
      % Cast to PlaneWave
      beam = castHelper(@ott.beam.vswf.PlaneWave.like, beam, varargin{:});
    end
  end

  methods (Sealed)
    function beam = ott.beam.abstract.Negative(beam, varargin)
      % Negate fields
      beam.field = -beam.field;
    end

    function beam = ott.beam.PlaneWave(beam, varargin)
      % Cast to PlaneWave
      beam = castArrayHelper(@ott.beam.PlaneWave.like, beam, varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Cast to Ray
      beam = castArrayHelper(@ott.beam.Ray.like, beam, varargin{:});
    end

    function beam = ott.beam.vswf.PlaneBasis(varargin)
      % Cast to PlaneWave
      %
      % Usage
      %   beam = ott.beam.vswf.PlaneBasis(beam, ...)
      %   Additional parameters are passed to vswf.PlaneBasis.like

      beam = castArrayHelper(@ott.beam.vswf.PlaneBasis.like, varargin{:});
    end
  end

  methods (Sealed, Access=protected)
    function beam = castArrayHelper(cast, beam, varargin)
      % Helper for casts

      assert(isa(beam, 'ott.beam.abstract.PlaneWaveStem'), ...
          'First argument must be a abstract.PlaneWaveStem');
      ott.utils.nargoutCheck(beam, nargout);

      % All casts use the PlaneWaveArray interface
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
    function beam = set.field(beam, val)
      assert(isnumeric(val) && ismatrix(val) && numel(val) == 2, ...
          'field must be 2 element numeric matrix');
      beam.field = val(:);
    end

    function beam = set.directionSet(beam, val)
      beam.rotation = val;
    end
    function dirset = get.directionSet(beam)
      dirset = beam.rotation;
    end

    function beam = set.origin(beam, val)
      beam.position = val;
    end
    function origin = get.origin(beam)
      origin = beam.position;
    end
  end
end

