classdef (Abstract) PlaneWaveArray < ott.beam.properties.PlaneWave ...
    & ott.beam.properties.ArrayType
% Properties for plane wave arrays.
% Inherits from :class:`PlaneWave` and :class:`ArrayType`.
%
% Properties
%   - field         -- Field parallel and perpendicular to polarisation
%   - directionSet  -- Set of direction vectors describing orientation
%   - origin        -- Position used to calculate beam phase offset
%   - position      -- (Inherited) Beam position
%   - rotation      -- (Inherited) Beam rotation
%
% Dependent properties
%   - direction     -- Beam direction
%   - polarisation1 -- Primary polarisation direction
%   - polarisation2 -- Secondary polarisation direction
%   - intensity     -- Intensity of the beam
%   - wavevector    -- (Inherited) Beam wave-vector
%
% Methods
%   - setData     -- Set field, directionSet and origin
%   - size        -- Size of PlaneWave array
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%   - translate*  -- Apply a phase offset (translation) to the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    field          % Field parallel and perpendicular to polarisation
    directionSet   % Set of direction vectors describing orientation
    origin         % Position used to calculate beam phase offset
  end

  properties (Dependent, SetAccess=protected)
    direction      % Beam direction
    polarisation1  % Primary polarisation direction
    polarisation2  % Secondary polarisation direction
    intensity      % Intensity of the beam
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.PlaneWaveArray')
        args = ott.utils.addDefaultParameter('field', other.field, args);
        args = ott.utils.addDefaultParameter('origin', other.origin, args);
        args = ott.utils.addDefaultParameter(...
            'directionSet', other.directionSet, args);
      end
      args = ott.beam.properties.PlaneWave.likeProperties(other, args);
    end

    function [origin, directionSet, field, unmatched] = parseArgs(varargin)
      % Helper for parsing arguments

      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('direction', [], @isnumeric);
      p.addOptional('polarisation1', [], @isnumeric);
      p.addParameter('directionSet', eye(3));
      p.addParameter('field', [1; 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      origin = p.Results.origin;
      field = p.Results.field;

      % Get directionSet from inputs
      hasDir = ~any(strcmpi('direction', p.UsingDefaults));
      hasPol1 = ~any(strcmpi('polarisation1', p.UsingDefaults));
      hasSet = ~any(strcmpi('directionSet', p.UsingDefaults));
      assert(hasSet + (hasDir & hasPol1) <= 1, ...
          'must only provide direction/polarisation or directionSet');
      if hasDir || hasPol1
        assert(hasDir && hasPol1, ...
            'must provide both direction and polarisation or directionSet');

        % Calculate secondary polarisation direction
        polarisation2 = cross(p.Results.direction, p.Results.polarisation1);

        % Form direction set
        directionSet = cat(2, ...
            reshape(p.Results.polarisation1, 3, 1, []), ...
            reshape(polarisation2, 3, 1, []), ...
            reshape(p.Results.direction, 3, 1, []));
        directionSet = reshape(directionSet, 3, []);
      else
        directionSet = p.Results.directionSet;
      end
    end
  end

  methods
    function bm = PlaneWaveArray(varargin)
      % Construct plane wave array properties
      %
      % Usage
      %   beam = PlaneWaveArray(origin, direction, polarisation1, ...)
      %
      %   beam = PlaneWaveArray(origin, 'directionSet', directionSet, ...)
      %
      % Parameters
      %   - origin (3xN numeric) -- Plane wave origins.
      %   - direction (3xN numeric) -- Plane wave direction.
      %   - polarisation1 (3xN numeric) -- Primary polarisation direction.
      %
      %   - directionSet (3x3N numeric) -- Array formed by combining
      %     direction/polarisation vectors into rotation matrices.  The
      %     direction vector should be the last column of the matrix.
      %     Default: ``eye(3)``.
      %
      % Optional named arguments
      %   - field (2xN numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %     Default: ``[1; 1i]``.

      [origin, directionSet, field, unmatched] = ...
          ott.beam.properties.PlaneWaveArray.parseArgs(varargin{:});

      bm = bm@ott.beam.properties.PlaneWave(unmatched{:});
      bm = bm.setData(origin, directionSet, field);
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

    function sz = size(beam, varargin)
      % Get the number of beams contained in this object
      %
      % Usage
      %   sz = size(beam)   or    sz = beam.size()
      %   For help on arguments, see builtin ``size``.
      %
      % The leading dimension is always 1.  May change in future.

      sz = size(beam.field);
      sz(1) = 1;

      [varargout{1:nargout}] = ott.utils.size_helper(sz, varargin{:});
    end
  end

  methods (Hidden)
    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams

      assert(dim == 2, 'Only 1xN arrays supported (may change in future)');

      other_origin = {};
      other_field = {};
      other_set = {};
      for ii = 1:length(varargin)
        other_origin{ii} = varargin{ii}.origin;
        other_field{ii} = varargin{ii}.field;
        other_set{ii} = varargin{ii}.directionSet;
      end

      beam.origin = cat(dim, beam.origin, other_origin{:});
      beam.field = cat(dim, beam.field, other_field{:});
      beam.directionSet = cat(dim, beam.directionSet, other_set{:});
    end

    function beam = plusInternal(beam1, beam2)
      % Concatenate two coherent beams together

      beam = cat(2, beam1, beam2);
    end

    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam

      if numel(subs) > ndims(beam.data)
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) > ndims(beam.origin), ...
            'Too many subscript indices');
      end

      % Must set data first!
      beam.origin = beam.origin(:, subs{:});
      beam.field = beam.field(:, subs{:});

      idx = (1:3).' + 3*(subs{1}-1);
      beam.directionSet = beam.directionSet(:, idx(:));
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam

      if numel(subs) > ndims(beam.data)
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) > ndims(beam.origin), ...
            'Too many subscript indices');
      end

      assert(isempty(rem), 'Assignment to parts of beams not supported');

      idx = (1:3).' + 3*(subs{1}-1);

      if isempty(other)
        % Delete data
        beam.origin(:, subs{:}) = other;
        beam.field(:, subs{:}) = other;
        beam.directionSet(:, idx(:)) = other;

      else
        % Ensure we have a plane wave
        % TODO: We used to have a cast here instead
        assert(isa(other, 'ott.beam.properties.PlaneWave'), ...
            'Only PlaneWave beams supported for now');

        % Must set data first!
        beam.origin(:, subs{:}) = other.origin;
        beam.field(:, subs{:}) = other.field;
        beam.directionSet(:, idx) = other.directionSet;
      end
    end
  end

  methods % Getters/setters
    % Properties
    %   - field         -- Field parallel and perpendicular to polarisation
    %   - directionSet  -- Set of direction vectors describing orientation
    %   - origin        -- Position used to calculate beam phase offset
    %   - direction     -- Beam direction
    %   - polarisation1 -- Primary polarisation direction
    %   - polarisation2 -- Secondary polarisation direction
    %   - intensity     -- Intensity of the beam

    function direction = get.direction(beam)
      direction = beam.directionSet(:, 3:3:end);
    end
    function direction = get.polarisation1(beam)
      direction = beam.directionSet(:, 1:3:end);
    end
    function direction = get.polarisation2(beam)
      direction = beam.directionSet(:, 2:3:end);
    end

    function intensity = get.intensity(beam)
      if strcmpi(beam.array_type, 'coherent')
        % TODO: Doesn't this need origin/rotations too?
        intensity = sum(abs(sum(beam.field, 2)), 1);
      elseif strcmpi(beam.array_type, 'incoherent')
        intensity = sum(sum(abs(beam.field)));
      else
        intensity = sum(abs(beam.field), 1);
      end
    end
  end
end
