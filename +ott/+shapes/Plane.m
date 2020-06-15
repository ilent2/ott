classdef Plane < ott.shapes.Shape ...
    & ott.shapes.mixin.CoordsCart ...
    & ott.shapes.mixin.InfVolume
% Shape describing a plane with infinite extent
% Inherits from :class:`ott.shapes.Shape`.
%
% Dependent properties
%   - normal      -- Vector representing surface normal
%   - offset      -- Offset of surface from coordinate origin
%
% Supported casts
%   - TriangularMesh    -- (Inherited) Uses PatchMesh
%   - PatchMesh
%   - Strata
%   - Slab

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    normal       % Vector representing surface normal
    offset       % Offset of surface from coordinate origin

    boundingBox        % Bounding box surrounding mesh (semi-infinite)
    starShaped         % Always false
    xySymmetry         % Always false
    zRotSymmetry       % Always 0
  end

  methods
    function shape = Plane(varargin)
      % Construct a new infinite plane
      %
      % Usage
      %   shape = Plane(normal, ...)
      %
      % Optional named arguments
      %   - normal (3xN numeric) -- Normals to planes.  Default: ``[]``.
      %     Overwrites any values set with `rotation`.
      %
      %   - offset (1xN numeric) -- Offset of the plane from the position.
      %     Default: ``[]``.  Overwrites any values set with `position`.
      %
      %   - position (3xN numeric) -- Position of the plane.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3N numeric) -- Plane orientations.
      %     Default: ``eye(3)``.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('normal', []);
      p.addParameter('offset', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shapes.Shape(unmatched{:});

      if ~isempty(p.Results.normal)
        shape.normal = p.Results.normal;
      end
      if ~isempty(p.Results.offset)
        shape.offset = p.Results.offset;
      end
    end

    function shape = ott.shapes.PatchMesh(shape, varargin)
      % Convert the plane into a patch
      %
      % Usage
      %   shape = ott.shapes.PatchMesh(shape, ...)
      %
      % Optional named parameters
      %   - scale (numeric) -- Size of the generated patch.
      %
      % Unmatched parameters are passed to patch constructor.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      v1 = shape.rotation(:, 1);
      v2 = shape.rotation(:, 2);

      v1 = v1 ./ vecnorm(v1) .* p.Results.scale ./ 2;
      v2 = v2 ./ vecnorm(v2) .* p.Results.scale ./ 2;

      verts = [v1 + v2, v1 - v2, -v1 - v2, -v1 + v2];
      faces = [1; 2; 3; 4];

      shape = ott.shapes.PatchMesh(verts, faces, ...
          'position', shape.position, 'rotation', eye(3), unmatched{:});
    end

    function shape = ott.shapes.Strata(planearray)
      % Array of planes can be cast to Strata if normals align

      % Check normals
      normals = [planearray.normal];
      assert(all(normals(1) == normals), 'all normals must match');

      % Calculate depth of each slab
      offsets = [shape.offset];
      depth = diff(offsets);

      % Create shape
      shape = ott.shapes.Strata(normal, depth, ...
          'offset', planearray(1).offset);
    end

    function shape = ott.shapes.Slab(planearray)
      % Array of two shapes can be cast to Slab

      stratashape = ott.shapes.Strata(planearray);
      shape = ott.shapes.Slab(stratashape);
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
      sz = size(vecs);
      norms = repmat(normal, [1, sz(2:end)]);

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

    function surf(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Converts the shape to a PatchMesh and calls surf.
      %
      % Usage
      %   shape.surf(...)
      %
      % Optional named parameters
      %   - scale (numeric) -- Size of patch.  Default: ``1.0``.
      %
      % Additional named parameters are passed to PatchMesh.surf.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shapes.PatchMesh(shape, 'scale', p.Results.scale);
      shape.surf(unmatched{:});
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz)
      % Determine if a point is on one side of the plane or the other

      % Determine if points are above plane
      b = (sum(xyz .* shape.normal, 1) - shape.offset) > 0;
    end

    function nxyz = normalsXyzInternal(shape, xyz)
      nxyz = repmat(shape.normal, 1, size(xyz, 2));
    end
  end

  methods % Getters/setters
    function plane = set.normal(plane, val)
      assert(isnumeric(val) && numel(val) == 3, ...
          'Normal must be 3 element numeric vector');

      val = val ./ vecnorm(val);

      n1 = cross(val, [0;0;1]);
      if vecnorm(n1) < 0.1
        n1 = cross(val, [1;0;0]);
      end
      n2 = cross(val, n1);

      plane.rotation = [n1./vecnorm(n1), n2./vecnorm(n2), val];
    end

    function n = get.normal(shape, val)
      n = shape.rotation(:, 3);
    end

    function plane = set.offset(plane, val)
      assert(isnumeric(val) && isscalar(val), ...
          'offset must be numeric scalar');
      plane.position = val .* plane.normal;
    end

    function o = get.offset(plane)
      o = dot(plane.position, plane.normal);
    end

    function bb = get.boundingBox(shape)
      bb = [-Inf, Inf; -Inf; Inf; -Inf; 0];
    end

    function b = get.starShaped(shape)
      b = false;
    end
    function b = get.xySymmetry(shape)
      b = false;
    end
    function q = get.zRotSymmetry(shape)
      q = 0;
    end
  end
end
