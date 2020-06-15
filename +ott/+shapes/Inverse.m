classdef Inverse < ott.shapes.Shape
% Inverts the geometry of a shape.
%
% Properties
%   - internal    -- Internal shape that is inversed
%   - volume      -- If volume was finite, makes it infinite
%   - maxRadius   -- If maxRadius was finite, makes it infinite
%
% Methods
%   - operator~   -- Smart inverse

  properties
    internal      % Internal shape that is unversed
  end

  properties (Dependent)
    volume
    maxRadius
    boundingBox
    starShaped        % Inherited from internal
    xySymmetry        % Inherited from internal
    zRotSymmetry      % Inherited from internal
  end

  methods
    function shape = Inverse(internal, varargin)
      % Construct a new inverse shape
      %
      % Usage
      %   shape = Inverse(internal, ...)
      %
      % Additional parameters passed to base.

      shape = shape@ott.shapes.Shape(varargin{:});
      shape.internal = internal;
    end

    function ishape = not(shape)
      % Take the inverse smartly

      % Apply rotation/translation to internal shape
      ishape = shape.internal.rotate(shape.rotation);
      ishape.position = ishape.position + shape.position;
    end

    function varargout = surf(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Usage
      %   p = shape.surf(...)
      %   Calls the corresponding `surf` method of the internal shape.

      [varargout{1:nargout}] = shape.internal.surf(varargin{:});
    end

    function [xyz, nxyz, dA] = surfPoints(shape, varargin)
      % Calculate points for surface integration.
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)
      %   Calls the corresponding method of the internal shape.

      [xyz, nxyz, dA] = shape.internal.surfPoints(varargin{:});
    end
  end

  methods (Hidden)
    function b = insideRtpInternal(shape, rtp)
      % Determine if point is inside inverse
      b = ~shape.internal.insideRtp(rtp, 'origin', 'global');
    end

    function b = insideXyzInternal(shape, xyz)
      % Determine if point is inside inverse
      b = ~shape.internal.insideXyz(xyz, 'origin', 'global');
    end

    function nxyz = normalsRtpInternal(shape, rtp)
      nxyz = -shape.internal.normalsRtp(rtp, 'origin', 'global');
    end

    function nxyz = normalsXyzInternal(shape, xyz)
      nxyz = -shape.internal.normalsXyz(xyz, 'origin', 'global');
    end
  end

  methods % Getters/setters
    function shape = set.internal(shape, val)
      assert(isa(val, 'ott.shapes.Shape') && numel(val) == 1, ...
          'internal must be a single ott.shapes.Shape');
      shape.internal = val;
    end

    function v = get.volume(shape)
      if isfinite(shape.internal.volume)
        v = Inf;
      else
        % TODO: numerical integration?
        warning('volume may not be infinite');
        v = Inf;
      end
    end

    function r = get.maxRadius(shape)
      if isfinite(shape.internal.maxRadius)
        r = Inf;
      else
        % TODO: numerical integration?
        warning('radius may not be infinite');
        r = Inf;
      end
    end

    function bb = get.boundingBox(shape)
      bb = [-Inf, Inf; -Inf, Inf; -Inf, Inf];
    end

    function b = get.starShaped(shape)
      b = shape.internal.starShaped;
    end
    function shape = set.starShaped(shape, val)
      shape.internal.starShaped = val;
    end

    function b = get.xySymmetry(shape)
      b = shape.internal.xySymmetry
    end
    function shape = set.xySymmetry(shape, val)
      shape.internal.xySymmetry = val;
    end

    function b = get.zRotSymmetry(shape)
      b = shape.internal.zRotSymmetry
    end
    function shape = set.zRotSymmetry(shape, val)
      shape.internal.zRotSymmetry = val;
    end
  end
end
