classdef RotationPositionProp < ott.utils.RotateHelper
% Adds rotation and position properties to a class
% Inherits from :class:`ott.utils.RotateHelper`.
%
% Properties
%   - position      -- 3x1 position of object
%   - rotation      -- 3x3 rotation of object

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    position = [0;0;0];
    rotation = eye(3);
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
