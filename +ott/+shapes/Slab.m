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
    function shape = Slab(normal, depth, varargin)
      % Construct a new infinite slab
      %
      % Usage
      %   shape = Slab(normal, depth, ...)
      %
      % Parameters
      %   - normal (3xN numeric) -- Surface normal
      %   - depth (N numeric) -- Depth of surface
      %
      % Optional named arguments
      %   - offset (N numeric) -- Offset of first surface to origin.
      %     Default: ``0``.
      %
      %   - position (3 numeric) -- Position of the shape.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Orientation of the shape.
      %     Default: ``eye(3)``.

      shape = shape@ott.shapes.Plane(normal, varargin{:});
      shape.depth = depth;
    end

    function shape = ott.shapes.Strata(oldshape)
      % Can be cast to a Strata

      shape = ott.shapes.Strata(oldshape.normal, oldshape.offset, ...
          oldshape.depth);
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
      shape.surfAddArgs(p);
      p.parse(varargin{:});

      % Calculate the X, Y, Z coordinates for a plane surface
      [X, Y, Z] = shape.calculateSurface(p);

      % Duplicate plane for second surface
      offset = shape.normal*(shape.offset+shape.depth);
      X = [X, nan(size(X, 1), 1), X + offset(1)];
      Y = [Y, nan(size(X, 1), 1), Y + offset(2)];
      Z = [Z, nan(size(X, 1), 1), Z + offset(3)];

      % Draw the figure and handle rotations/translations
      [varargout{1:nargout}] = shape.surfCommon(p, size(X), X, Y, Z);
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz)
      % Determine if a point is on one side of the plane or the other

      % Determine if points are inside slab
      dist = sum(xyz .* shape.normal, 1) - shape.offset;
      b = dist > 0 & (dist - shape.depth < 0);
    end
  end
end
