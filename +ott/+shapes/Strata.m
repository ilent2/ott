classdef Strata < ott.shapes.Plane
% Shape describing a series of stratified interfaces
%
% This shape describes a series of layered planes.  When the number
% of layers is equal to 2, this object can be converted to a Slab.
% All points above the first layer are considered to be inside the shape.
%
% Properties
%   - normal      -- Vector representing surface normal
%   - offset      -- Offset of surface from coordinate origin
%   - depth       -- Depth of each layer

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    depth         % Depth of each layer
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
      %   - depth (numeric) -- Depth of each layer.
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

    function shape = ott.shapes.Slab(oldshape)
      % Can be cast to a slab if the number of surfaces is 2

      assert(numel(oldshape.depth) == 1, 'Depth must have 1 element');

      shape = ott.shapes.Slab(oldshape.normal, ...
          oldshape.offset, oldshape.depth);
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

      % Calculate the X, Y, Z coordinates for a plane surface
      [X0, Y0, Z0] = shape.calculateSurface(p)

      % Add plane for first interface
      X = X0;
      Y = Y0;
      Z = Z0;

      % Duplicate plane for remaining interfaces
      for ii = 1:length(depth)
        offset = shape.normal*(shape.offset+sum(shape.depth(1:ii)));
        X = [X, nan(size(X, 1), 1), X0 + offset(1)];
        Y = [Y, nan(size(X, 1), 1), Y0 + offset(2)];
        Z = [Z, nan(size(X, 1), 1), Z0 + offset(3)];
      end

      % Draw the figure and handle rotations/translations
      shape.surfCommon(p, sz, xx, yy, zz);
    end
  end
end
