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
%
% Casts
%   - ott.beam.vswf.PlaneBasis
%   - Ray
%   - PlaneWave

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    field          % Field parallel and perpendicular to polarisation
    directionSet   % Set of direction vectors describing orientation
    origin         % Position used to calculate beam phase offset
  end

  methods
    function varargout = size(beam, varargin)
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

      if numel(subs) > ndims(beam.origin)
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) > ndims(beam.origin), ...
            'Too many subscript indices');
      end

      beam.origin = beam.origin(:, subs{:});
      beam.field = beam.field(:, subs{:});

      idx = (1:3).' + 3*(subs{1}-1);
      beam.directionSet = beam.directionSet(:, idx(:));
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam

      if numel(subs) > ndims(beam.origin)
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
end
