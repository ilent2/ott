classdef Slab < ott.shapes.Plane
% Shape describing a slab with infinite extent in two directions
% Inherits from :class:`ott.shapes.Plane`.
%
% Properties
%   - normal      -- Vector representing surface normal
%   - offset      -- Offset of surface from coordinate origin
%   - depth       -- Depth of the slab

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    depth         % Depth of the slab
  end

  methods
    function shape = Slab(normal, offset, depth)
      % Construct a new infinite slab
      %
      % Usage
      %   shape = Slab(normal, depth)
      %
      %   shape = Slab(normal, offset, depth)
      %
      % Parameters
      %   - normal (3x1 numeric) -- Surface normal
      %   - depth (numeric) -- Depth of surface
      %   - offset (numeric) -- Offset of first surface to origin.

      if nargin == 2
        shapeArgs = {normal};
        depth = offset;
      elseif nargin == 3
        shapeArgs = {normal, offset};
      else
        error('Must supply 2 or 3 inputs arguments');
      end

      shape = shape@ott.shapes.Plane(shapeArgs{:});
      shape.depth = depth;
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

      % Determine if points are inside slab
      dist = sum(xyz .* shape.normal, 1) - shape.offset;
      b = dist > 0 & (dist - shape.depth < 0);
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
      data = [data, nan(3, 2), [v1, v2, v3, v4] ...
          + shape.normal*(shape.offset+shape.depth)];
      X = reshape(data(1, :), [2, 5]);
      Y = reshape(data(2, :), [2, 5]);
      Z = reshape(data(3, :), [2, 5]);

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
