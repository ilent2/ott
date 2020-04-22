classdef Plane < ott.shapes.ShapeCart
% Shape describing a plane with infinite extent
% Inherits from :class:`ott.shapes.ShapeCart`.
%
% Properties
%   - normal      -- Vector representing surface normal
%   - offset      -- Offset of surface from coordinate origin

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    normal       % Vector representing surface normal
    offset       % Offset of surface from coordinate origin
  end

  methods
    function shape = Plane(normal, offset)
      % Construct a new infinite plane
      %
      % Usage
      %   shape = Plane(normal)
      %
      %   shape = Plane(normal, offset)

      shape = shape@ott.shapes.ShapeCart();

      if nargin == 1
        offset = 0.0;
      end

      shape.normal = normal;
      shape.offset = offset;
    end

    function r = get_maxRadius(shape)
      % Infinite plane has infinite maximum radius
      r = Inf;
    end

    function v = get_volume(shape)
      % Infinite plane has infinite volume
      v = Inf;
    end

    function b = insideXyz(shape, varargin)
      % Determine if a point is on one side of the plane or the other
      %
      % Usage
      %   b = shape.insideXyz(x, y, z, ...) determine if the point
      %   is above or bellow the plane surface.
      %
      %   b = shape.insideXyz(xyz, ...) as above but with a 3xN matrix.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.

      % Get xyz coordinates from inputs and translated to origin
      xyz = shape.insideXyzParseArgs(shape.position, varargin{:});

      % Determine if points are above plane
      b = (sum(xyz .* shape.normal, 1) - shape.offset) > 0;
    end

    function varargout = surf(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Usage
      %   shape.surf(...) displays a visualisation of the shape in
      %   the current figure.
      %
      %   [X, Y, Z] = shape.surf() calculates the coordinates and
      %   arranges them in a grid suitable for use with matlab surf function.
      %
      % Optional named arguments
      %   - scale (numeric) -- Scaling factor for the plane.
      %
      %   - axes ([] | axes handle) -- axis to draw in.  Default: ``gca``.
      %
      %   - surfoptions (cell array) -- options to be passed to surf.
      %     Default: ``{}``.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.addParameter('surfoptions', {});
      p.addParameter('axes', []);
      p.parse(varargin{:});

      % Calculate two orthogonal vectors
      v = shape.normal;
      [~, I] = min(abs(v));
      v(I) = max(abs(v));
      c1 = cross(v, shape.normal);
      c2 = cross(c1, shape.normal);

      % Apply scale to c1 and c2
      c1 = c1 ./ vecnorm(c1) .* p.Results.scale;
      c2 = c2 ./ vecnorm(c2) .* p.Results.scale;

      % Calculate 4 corners of plane
      v1 = c1 + c2 + shape.normal*shape.offset;
      v2 = c1 - c2 + shape.normal*shape.offset;
      v3 = -c1 + c2 + shape.normal*shape.offset;
      v4 = -c1 - c2 + shape.normal*shape.offset;

      % Reshape into surf form
      data = [v1, v2, v3, v4];
      X = reshape(data(1, :), [2, 2]);
      Y = reshape(data(2, :), [2, 2]);
      Z = reshape(data(3, :), [2, 2]);

      % Generate the surface
      if nargout == 0 || ~isempty(p.Results.axes)

        % Place the surface in the specified axes
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        surf(our_axes, X, Y, Z, p.Results.surfoptions{:});
      end

      % Set outputs if requested
      if nargout ~= 0
        varargout = { X, Y, Z };
      end
    end
  end
end
