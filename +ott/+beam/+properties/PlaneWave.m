classdef (Abstract) PlaneWave < ott.beam.properties.Beam ...
    & ott.beam.properties.VariableMedium
% Properties of scattered beams.
% Inherits from :class:`Beam`.
%
% Static methods
%   - likeProperties    -- Construct argument list of like-properties
%   - DirectionSet      -- Construct directionSet
%
% Properties
%   - position      -- (Inherited) Beam position
%   - rotation      -- (Inherited) Beam rotation
%
% Abstract properties
%   - field         -- Field parallel and perpendicular to polarisation
%   - origin        -- Position used to calculate beam phase offset
%   - directionSet  -- Set of directions describing 
%
% Dependent properties
%   - wavevector    -- Beam wave-vector
%   - direction     -- Beam direction
%   - polarisation1 -- Primary polarisation direction
%   - polarisation2 -- Secondary polarisation direction
%   - intensity     -- Intensity of the beam
%
% Methods
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%   - translate*  -- Apply a phase offset (translation) to the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract, SetAccess=protected)
    field             % Field parallel and perpendicular to polarisation
    origin            % Position used to calculate beam phase offset
    directionSet      % Set of direction vectors describing orientation
  end

  properties (Dependent)
    wavevector        % Wave-vectors of plane wave components
    direction         % Beam direction
    polarisation1     % Primary polarisation direction
    polarisation2     % Secondary polarisation direction
    intensity         % Intensity of plane wave components
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.PlaneWave')
        args = ott.utils.addDefaultParameter('field', other.field, args);
        args = ott.utils.addDefaultParameter('origin', other.origin, args);
        args = ott.utils.addDefaultParameter(...
            'directionSet', other.directionSet, args);
      end
      args = ott.beam.properties.Beam.likeProperties(other, args);
      args = ott.beam.properties.VariableMedium.likeProperties(other, args);
    end

    function directionSet = DirectionSet(direction, polarisation1)
      % Construct a new direction set from direction/polarisation
      %
      % Usage
      %   directionSet = DirectionSet(direction, polarisation1)
      %
      % Parameters
      %   - direction (3xN numeric) -- Propagation direction of wave.
      %   - polarisation (3xN numeric) -- Primary polarisation direction.

      assert(isnumeric(direction) && ismatrix(direction) ...
          && size(direction, 1) == 3, ...
          'direction must be 3xN numeric matrix');
      assert(isnumeric(polarisation1) && ismatrix(polarisation1) ...
          && all(size(polarisation1) == size(direction)), ...
          'polarisation1 must be 3xN numeric with same size as direction');
        
      direction = direction ./ vecnorm(direction, 2, 1);
      polarisation1 = polarisation1 ./ vecnorm(polarisation1, 2, 1);

      % Calculate secondary polarisation direction
      polarisation2 = cross(direction, polarisation1);

      % Form direction set
      directionSet = cat(2, ...
          reshape(polarisation1, 3, 1, []), ...
          reshape(polarisation2, 3, 1, []), ...
          reshape(direction, 3, 1, []));
      directionSet = reshape(directionSet, 3, []);
    end
  end

  methods
    function beam = PlaneWave(varargin)
      % Construct a new plane wave properties description
      %
      % Usage
      %   beam = PlaneWave(origin, directionSet, field, ...)
      %   Parameters can also be named arguments.
      %
      % See also :meth:`FromDirection` for alternative constructor
      % using direction/polarisation.
      %
      % Parameters
      %   - origin (3xN numeric) -- Plane wave origins.
      %
      %   - directionSet (3x3N numeric) -- Array formed by combining
      %     direction/polarisation vectors into rotation matrices.  The
      %     direction vector should be the last column of the matrix.
      %
      %   - field (2xN numeric) -- Field in two orthogonal polarisation
      %     directions (described by first two columns of directionSet).

      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('directionSet', [], @isnumeric);
      p.addOptional('field', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Beam(unmatched{:});
      beam = beam.setData(p.Results.origin, ...
          p.Results.directionSet, p.Results.field);
    end

    function beam = setData(beam, origin, directionSet, field)
      % Set beam data
      %
      % Usage
      %   beam = beam.setData(origin, directionSet, field)
      %
      % Parameters
      %   - origin (3xN numeric) -- Plane wave origins.
      %
      %   - directionSet (3x3N numeric) -- Array formed by combining
      %     direction/polarisation vectors into rotation matrices.  The
      %     direction vector should be the last column of the matrix.
      %     Default: ``eye(3)``.
      %
      %   - field (2xN numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %
      % Check implementation for supported array sizes.

      ott.utils.nargoutCheck(beam, nargout);

      Norigin = size(origin, 2);
      Nset = size(directionSet, 2)./3;
      Nfield = size(field, 2);
      Nelem = max([Norigin, Nset, Nfield]);

      assert(isnumeric(field) && ismatrix(field) && size(field, 1) == 2, ...
          'field must be 2xN numeric matrix');
      assert(isnumeric(origin) && ismatrix(origin) && size(origin, 1) == 3, ...
          'origin must be 3xN numeric matrix');
      assert(isnumeric(directionSet) && ismatrix(directionSet) ...
          && size(directionSet, 1) == 3 ...
          && mod(size(directionSet, 2), 3) == 0, ...
          'directionSet must be 3x3N numeric matrix');
      assert(Norigin == Nelem || Norigin == 1, ...
          'origin must be length 1 or same length as other arguments');
      assert(Nset == Nelem || Nset == 1, ...
          'directionSet must be length 1 or same as other inputs');
      assert(Nfield == Nelem || Nfield == 1, ...
          'field must be length 1 or same as other inputs');

      % Duplicate scalars
      if Norigin == 1, origin = repmat(origin, 1, Nelem); end
      if Nfield == 1, field = repmat(field, 1, Nelem); end
      if Nset == 1, directionSet = repmat(directionSet, 1, Nelem); end

      beam.origin = origin;
      beam.directionSet = directionSet;
      beam.field = field;
    end
  end

  methods % Getters/setters
    % Properties
    %   wavevector        % Wave-vectors of plane wave components
    %   direction         % Beam direction
    %   polarisation1     % Primary polarisation direction
    %   polarisation2     % Secondary polarisation direction
    %   intensity         % Intensity of plane wave components

    function wv = get.wavevector(beam)
      % Get the plane wave wave-vector
      wv = beam.direction .* beam.wavenumber ./ vecnorm(beam.direction);
    end

    function direction = get.direction(beam)
      direction = beam.directionSet(:, 3:3:end);
    end
    function beam = set.direction(beam, val)
      assert(~isempty(val), 'direction must not be empty');
      beam.directionSet(:, 3:3:end) = val;
    end

    function direction = get.polarisation1(beam)
      direction = beam.directionSet(:, 1:3:end);
    end
    function beam = set.polarisation1(beam, val)
      assert(~isempty(val), 'polarisation1 must not be empty');
      beam.directionSet(:, 1:3:end) = val;
    end

    function direction = get.polarisation2(beam)
      direction = beam.directionSet(:, 2:3:end);
    end
    function beam = set.polarisation2(beam, val)
      assert(~isempty(val), 'polarisation2 must not be empty');
      beam.directionSet(:, 2:3:end) = val;
    end

    function intensity = get.intensity(beam)
      intensity = sum(abs(beam.field), 1);
    end
    function beam = set.intensity(beam, val)
      beam.field = beam.field ./ beam.intensity .* val;
    end
  end
end
