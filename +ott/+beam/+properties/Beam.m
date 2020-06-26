classdef (Abstract) Beam < ott.utils.RotationPositionProp
% A base class for Beam and abstract.Beam representations.
% Inherits from :class:`ott.utils.RotationPositionProp`.
%
% Any units can be used for the properties as long as they are
% consistent in all specified properties.  Calculated quantities
% will have these units.
%
% This class defines the common properties and methods to these
% two classes.
%
% Properties
%   - position      -- Position of the beam or array
%   - rotation      -- Rotation of the beam or array
%
% Abstract properties
%   - power         -- The power of the beam (may be infinite)
%   - medium        -- Medium where beam is propagating
%
% Methods
%  - rotate*     -- Beam rotation methods
%  - translate*  -- Beam translation methods

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract)
    power           % The power of the beam (may be infinite)
    medium          % Medium where beam is propagating
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      %
      % Usage
      %   args = Beam.like(other, args)

      args = ott.utils.addDefaultParameter('position', other.position, args);
      args = ott.utils.addDefaultParameter('rotation', other.rotation, args);
      args = ott.utils.addDefaultParameter('medium', other.medium, args);
    end
  end

  methods
    function beam = Beam(varargin)
      % Initialize beam properties
      %
      % Usage
      %   beam = beam@ott.beam.Properties(...)
      %
      % Named arguments
      %   - medium (ott.beam.medium.Medium|Material) -- Medium or
      %     material describing optical properties of medium.
      %     Default: ``ott.beam.medium.Vacuum.Unitary``.
      %
      %   - power (numeric) -- Initial beam power (if supported).
      %     No default, only sets property if argument passed.
      %
      %   - position (3 numeric) -- Position of beam.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Rotation matrix for beam.
      %     Default: ``eye(3)``.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.addParameter('power', []);
      p.addParameter('medium', ott.beam.medium.Vacuum.Unitary);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Store power if supplied
      if ~any(strcmpi('power', p.UsingDefaults))
        beam.power = p.Results.power;
      end

      % Store medium if supplied
      if ~any(strcmpi('medium', p.UsingDefaults))
        beam.medium = p.Results.medium;
      end

      % Store position/rotation
      beam.position = p.Results.position;
      beam.rotation = p.Results.rotation;
    end

    function b = contains(beam, array_type)
      % Query if a array_type is contained in the array.
      %
      % Usage
      %   b = beam.contains(array_type)
      %
      % Parameters
      %   - array_type (enum) -- An array type, must be one of
      %     'array', 'coherent' or 'incoherent'.

      % Terminal case: not an array.
      b = false;
    end

    function beam = defaultArrayType(beam, array_type, elements)
      % Construct a new array for this type

      if beam.contains('array')
        assert(strcmpi(array_type, 'array'), ...
            'type must be array for beams containing generic arrays');
      elseif beam.contains('incoherent')
        assert(any(strcmpi(array_type, {'array', 'incoherent'})), ...
            'type must be array/incoherent for incoherent content');
      end

      if nargin == 2
        beam = @(arg) ott.beam.Array(array_type, arg);
      else
        beam = ott.beam.Array(array_type, elements);
      end
    end

    function data = arrayApply(beam, func, varargin)
      % Apply function to each array in the beam array output.
      %
      % Usage
      %   data = beam.arrayApply(func, ...)
      %   Additional parameters are passed to the function.
      %
      % This function is overloaded by Array types in order to
      % implement incoherent combination.

      data = func(varargin{:});
    end
  end
end
