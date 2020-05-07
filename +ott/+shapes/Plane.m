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
    function shape = Plane(normal, varargin)
      % Construct a new infinite plane
      %
      % Usage
      %   shape = Plane(normal, ...)
      %
      % Parameters
      %   - normal (3xN numeric) -- Normals to planes.
      %
      % Optional named arguments
      %   - position (3xN numeric) -- Position of the plane.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3N numeric) -- Plane orientations.
      %     Default: ``eye(3)``.
      %
      %   - offset (1xN numeric) -- Offset of the plane from the position.
      %     Default: ``0.0``.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('offset', 0.0);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shapes.ShapeCart(unmatched{:});

      Nnormal = size(normal, 2);
      Noffset = numel(p.Results.offset);

      assert(Nnormal == 1 || numel(shape) == 1 || Nnormal == numel(shape), ...
          'length of normal must match length of position');
      assert(Noffset == 1 || numel(shape) == 1 || Noffset == numel(shape), ...
          'length of offset must match length of position');

      if numel(shape) == 1 && (Nnormal ~= 1 || Noffset ~= 1)
        shape = repelem(shape, 1, max([Nnormal, Noffset]));
      end

      % Assign normal
      if Nnormal > 1
        for ii = 1:Nnormal
          shape(ii).normal = normal(:, ii);
        end
      else
        [shape.normal] = deal(normal);
      end

      % Assign offset
      if Noffset > 1
        for ii = 1:Noffset
          shape(ii).offset = p.Results.offset(ii);
        end
      else
        [shape.offset] = deal(p.Results.offset);
      end
    end

    function shape = ott.shapes.Strata(planearray)
      % Array of planes can be cast to Strata if normals align

      % Check normals
      normal = planearray(1).normal;
      for ii = 2:numel(planearray)
        assert(all(normal == planearray(ii).normal), ...
            'All normals must match');
      end

      % Calculate depth of each slab
      offsets = [shape.offset];
      depth = diff(offsets);

      % Create shape
      shape = ott.shapes.Strata(normal, planearray(1).offset, depth);
    end

    function shape = ott.shapes.Slab(planearray)
      % Array of two shapes can be cast to Slab

      stratashape = ott.shapes.Strata(planearray);
      shape = ott.shapes.Slab(stratashape);
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

    function [locs, norms] = intersect(shape, vecs)
      % Calculate the intersection point on the plane surface.
      %
      % Rays will intersect the plane as long as they are traveling
      % towards the surface.  The normal will always be the surface
      % normal.
      %
      % Usage
      %   [locs, norms] = shape.intersect(vec)
      %   Returns a 3xN matrix of intersection locations or nan.
      %
      % Parameters
      %   - vec (utils.Vector) -- A vector or type that can be cast
      %     to a Vector.

      % Apply rotation to normal
      normal = shape.rotation * shape.normal;

      % Duplicate the normals
      norms = repmat(normal, 1, numel(vecs));

      % Calculate intersection location relative to ray origin
      dirs = vecs.direction ./ vecnorm(vecs.direction);
      ndirs = dot(dirs, norms);
      locs = -(dot(vecs.origin - shape.position, norms) ...
          - shape.offset) .* dirs ./ ndirs;

      % Remove rays traveling away from the plane
      locs(:, dot(locs, dirs) < 0) = nan;

      % Translate to vector origin
      locs = vecs.origin + locs;

      % Ensure infs are also nans
      locs(~isfinite(locs)) = nan;
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

      % Draw the figure and handle rotations/translations
      [varargout{1:nargout}] = shape.surfCommon(p, size(X), X, Y, Z);
    end
  end

  methods (Hidden)
    function surfAddArgs(beam, p)
      % Add surface drawing args to the input parser for surf
      p.addParameter('scale', 1.0);
      surfAddArgs@ott.shapes.ShapeCart(beam, p);
    end
    
    function [X, Y, Z] = calculateSurface(shape, p)
      % Calculate the X, Y, Z coordinates for a plane surface

      % Apply rotation to normal
      normal = shape.rotation * shape.normal;

      % Calculate two orthogonal vectors
      v = normal;
      [~, I] = min(abs(v));
      v(I) = max(abs(v));
      c1 = cross(v, normal);
      c2 = cross(c1, normal);

      % Apply scale to c1 and c2
      c1 = c1 ./ vecnorm(c1) .* p.Results.scale;
      c2 = c2 ./ vecnorm(c2) .* p.Results.scale;

      % Calculate 4 corners of plane
      v1 = c1 + c2 + normal*shape.offset;
      v2 = c1 - c2 + normal*shape.offset;
      v3 = -c1 + c2 + normal*shape.offset;
      v4 = -c1 - c2 + normal*shape.offset;

      % Reshape into surf form
      data = [v1, v2, v3, v4];
      X = reshape(data(1, :), [2, 2]) + shape.position(1);
      Y = reshape(data(2, :), [2, 2]) + shape.position(2);
      Z = reshape(data(3, :), [2, 2]) + shape.position(3);
    end
  end

  methods % Getters/setters
    function plane = set.normal(plane, val)
      assert(isnumeric(val) && numel(val) == 3, ...
          'Normal must be 3 element numeric vector');
      plane.normal = val(:);
    end

    function plane = set.offset(plane, val)
      assert(isnumeric(val) && isscalar(val), ...
          'offset must be numeric scalar');
      plane.offset = val;
    end
  end
end
