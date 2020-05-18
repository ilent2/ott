classdef (Abstract) RotateHelper
% A helper class providing may short-hand rotation methods.
%
% This class ensures a consistent rotation interface for all classes
% implementing rotations in the toolbox.  To use the class, simply
% inherit from this class and implement the ``rotate`` method.
%
% The toolbox uses the convention that angles are applied in the
% order X, Y, Z.  Most functions only use Y, Z.  This class only
% implements methods in this order.
%
% Methods
%   - rotate      -- Rotate by 3x3 rotation matrix
%   - rotateX     -- Rotate about x
%   - rotateY     -- Rotate about y
%   - rotateZ     -- Rotate about z
%   - rotateXy    -- Rotate about x then y
%   - rotateXz    -- rotate about x then z
%   - rotateYz    -- Rotate about y then z
%   - rotateXyz   -- Rotate about x, y and then z
%
% Abstract methods
%   - rotateInternal    -- The rotation method called by all other methods

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Abstract)
    rotateInternal     % The rotation method called by all other methods
  end

  methods
    function varargout = rotate(obj, mat, varargin)
      % Rotate by a 3x3 rotation matrix
      %
      % Usage
      %   [obj, ...] = obj.rotate(mat, ...)
      %   Applies the 3x3 rotation matrix.

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotateInternal(mat);
    end

    function varargout = rotateX(obj, angle, varargin)
      % Rotate the obj about the x axis
      %
      % Usage
      %   [obj, ...] = obj.rorateX(anglex, ...)
      %   Angles are specified in radians.

      import ott.utils.rotx;

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotate(rotx(angle*180/pi), varargin{:});
    end

    function varargout = rotateY(obj, angle, varargin)
      % Rotate the obj about the y axis
      %
      % Usage
      %   [obj, ...] = obj.rorateY(angley, ...)
      %   Angles are specified in radians.

      import ott.utils.roty;

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotate(roty(angle*180/pi), varargin{:});
    end

    function varargout = rotateZ(obj, angle, varargin)
      % Rotate the obj about the z axis
      %
      % Usage
      %   [obj, ...] = obj.rorateZ(anglez, ...)
      %   Angles are specified in radians.

      import ott.utils.rotz;

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotate(rotz(angle*180/pi), varargin{:});
    end

    function varargout = rotateXy(obj, anglex, angley, varargin)
      % Rotate the obj about the x then y
      %
      % Usage
      %   [obj, ...] = obj.rorateXy(anglex, angley, ...)
      %   Angles are specified in radians.

      import ott.utils.rotx;
      import ott.utils.roty;

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotate(roty(angley*180/pi) ...
          *rotx(anglex*180/pi), varargin{:});
    end

    function varargout = rotateXz(obj, anglex, anglez, varargin)
      % Rotate the obj about the x then z axes
      %
      % Usage
      %   [obj, ...] = obj.rorateXz(anglex, anglez, ...)
      %   Angles are specified in radians.
      %   Additional arguments are forwarded to rotate.

      import ott.utils.rotx;
      import ott.utils.rotz;

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotate(rotz(anglez*180/pi) ...
          *rotx(anglex*180/pi), varargin{:});

      obj.nargoutCheck(nargout);
    end

    function varargout = rotateYz(obj, angley, anglez, varargin)
      % Rotate the obj about the y then z axes
      %
      % Usage
      %   [obj, ...] = obj.rorateYz(angley, anglez, ...)
      %   Angles are specified in radians.
      %   Additional arguments are forwarded to rotate.

      import ott.utils.roty;
      import ott.utils.rotz;

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotate(rotz(anglez*180/pi) ...
          *roty(angley*180/pi), varargin{:});
    end

    function varargout = rotateXyz(obj, anglex, angley, anglez, varargin)
      % Rotate the obj about the x, y then z axes
      %
      % Usage
      %   [obj, ...] = obj.rorateXyz(anglex, angley, anglez, ...)
      %   Angles are specified in radians.
      %   Additional arguments are forwarded to rotate.

      import ott.utils.rotx;
      import ott.utils.roty;
      import ott.utils.rotz;

      ott.utils.nargoutCheck(obj, nargout);

      [varargout{1:nargout}] = obj.rotate(rotz(anglez*180/pi)* ...
          roty(angley*180/pi)*rotx(anglex*180/pi), varargin{:});
      obj.nargoutCheck(nargout);
    end
  end
end

