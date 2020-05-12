classdef RotationPositionProp < ott.utils.RotateHelper
% Adds rotation and position properties to a class
% Inherits from :class:`ott.utils.RotateHelper`.
%
% Properties
%   - position      -- 3x1 position of object
%   - rotation      -- 3x3 rotation of object
%
% Methods
%   - rotate*           -- Apply rotations to object
%   - forwardTransform  -- Coordinate transform: local to global
%   - inverseTransform  -- Coordinate transform: global to local
%
% Static methods
%   - rotationPositionHelper -- Helper for processing object arrays

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    position = [0;0;0];
    rotation = eye(3);
  end

  methods (Static)
    function varargout = rotationPositionHelper(func, other, varargin)
      % Helper function for processing arrays of objects.
      %
      % Dispatches calls to a function after translating or rotating
      % a particle by a specified position/rotation or an array of
      % positions and rotations.
      %
      % Usage
      %   [varargout{1:nargout}] = rotationPositionHelper(func, other, ...)
      %
      % Parameters
      %   - func (function_handle) -- Function to call with `other` and
      %     any unmatched arguments.
      %
      %   - other -- An object with rotation and position properties.
      %
      % Named arguments
      %   - position (3xN numeric) -- Position to apply to other before
      %     dispatching to func.  Default: ``[]``.
      %
      %   - rotation (3x3N numeric) -- Rotation to apply to other before
      %     dispatching to func.  Default: ``[]``.
      %
      %   - prxcat (numeric|empty) -- Dimension to concatenate prxfun
      %     results along.  Default: ``[]``.
      %
      % Length of position and rotation must match or be equal to 1.

      % TODO: Should this be its own function or merge with prxfun?

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('prxcat', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Check if we need to operate on multiple beams
      if size(p.Results.position, 2) > 1 || size(p.Results.rotation, 2) > 3
        [varargout{1:nargout}] = ott.utils.prxfun(@(varargin) ...
            ott.utils.RotationPositionProp.rotationPositionHelper(...
            func, other, varargin{:}), ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, unmatched{:});

        % Concatenate outputs along requested dimension
        if ~isempty(p.Results.prxcat)
          for ii = 1:nargout
            varargout{ii} = cat(p.Results.prxcat, varargout{ii}{:});
          end
        end

        return;
      end

      % Add position/rotation to beam
      if ~isempty(p.Results.position)
        other.position = other.position - p.Results.position;
      end
      if ~isempty(p.Results.rotation)
        other.rotation = p.Results.rotation.' * other.rotation;
      end

      [varargout{1:nargout}] = func(other, unmatched{:});
    end
  end

  methods
    function shape = rotate(shape, R)
      % Apply the rotation matrix to the shapes internal rotation

      % Check at least one output
      shape.nargoutCheck(nargout);

      % Apply rotation
      shape = R * shape.rotation;
    end

    function [xyz, vxyz] = forwardTransform(obj, xyz, vxyz)
      % Apply the position and rotation transformations to the coordinates
      %
      % Usage
      %   xyz = obj.forwardTransform(xyz)
      %   Apply transformation to position coordinates.
      %
      %   [xyz, vxyz] = obj.applyLocality(xyz, vxyz)
      %   Applies transformations to the position and vector coordinates.
      %
      % The `rotation` is applied to both `xyz` and `vxyz`.
      % The `position` is only applied to `xyz`.

      xyz = obj.rotation * xyz;
      xyz = xyz + obj.position;
      varargout{1} = xyz;

      if nargin == 3
        assert(nargout == 2, 'Too few output arguments');
        varargout{2} = obj.rotation * vxyz;
      else
        assert(nargout == 1, 'Too many output arguments');
      end
    end

    function varargout = inverseTransform(obj, xyz, vxyz)
      % Apply the position and rotation transformations to the coordinates
      %
      % Usage
      %   xyz = obj.forwardTransform(xyz)
      %   Apply transformation to position coordinates.
      %
      %   [xyz, vxyz] = obj.applyLocality(xyz, vxyz)
      %   Applies transformations to the position and vector coordinates.
      %
      % The `rotation` is applied to both `xyz` and `vxyz`.
      % The `position` is only applied to `xyz`.

      xyz = xyz - obj.position;
      xyz = obj.rotation.' * xyz;
      varargout{1} = xyz;

      if nargin == 3
        assert(nargout == 2, 'Too few output arguments');
        varargout{2} = obj.rotation.' * vxyz;
      else
        assert(nargout == 1, 'Too many output arguments');
      end
    end
  end

  methods % Getters/setters
    function obj = set.position(obj, val)
      assert(isnumeric(val) && all(size(val) == [3, 1]), ...
          'position must be a 3x1 numeric vector');
      obj.position = val;
    end

    function obj = set.rotation(obj, val)
      assert(isnumeric(val) && all(size(val) == [3, 3]), ...
          'rotation must be a 3x3 numeric vector');
      obj.rotation = val;
    end
  end
end
