classdef Slab < ott.shapes.Shape ...
    & ott.shapes.mixin.CoordsCart ...
    & ott.shapes.mixin.InfVolume
% Shape describing a slab with infinite extent in two directions
% Inherits from :class:`ott.shapes.Shape`.
%
% Properties
%   - normal      -- Vector representing surface normal
%   - depth       -- Depth of the slab
%
% Supported casts
%   - TriangularMesh    -- (Inherited) Uses PatchMesh
%   - PatchMesh         -- Uses Strata
%   - Strata

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    depth         % Depth of the slab
  end

  properties (Dependent)
    normal             % Vector representing surface normal
    boundingBox        % Bounding box surrounding mesh (semi-infinite)
    starShaped         % Always false
    xySymmetry         % Always false
    zRotSymmetry       % Always 0
  end

  methods
    function shape = Slab(varargin)
      % Construct a new infinite slab
      %
      % Usage
      %   shape = Slab(normal, depth, ...)
      %
      % Optional named arguments
      %   - depth (N numeric) -- Depth of surface.  Default: 0.5.
      %
      %   - normal (3xN numeric) -- Normals to planes.  Default: ``[]``.
      %     Overwrites any values set with `rotation`.
      %
      %   - position (3xN numeric) -- Position of the plane.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3N numeric) -- Plane orientations.
      %     Default: ``eye(3)``.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('normal', []);
      p.addOptional('depth', 0.5);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shapes.Shape(unmatched{:});
      shape.depth = p.Results.depth;

      if ~isempty(p.Results.normal)
        shape.normal = p.Results.normal;
      end
    end

    function shape = ott.shapes.PatchMesh(shape, varargin)
      % Cast to PatchMesh via Strata
      %
      % Usage
      %   shape = ott.shapes.PatchMesh(shape, ...)
      %
      % Optional named parameters
      %   - scale (numeric) -- Scale of the patch. Default: 1.0.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shapes.Strata(shape, unmatched{:});
      shape = ott.shapes.PatchMesh(shape, 'scale', p.Results.scale);
    end

    function shape = ott.shapes.Strata(shape, varargin)
      % Can be cast to a Strata

      shape = ott.shapes.Strata('normal', shape.normal, ...
          'depths', shape.depth, ...
          'position', shape.position-shape.depth./2, varargin{:});
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz)
      % Determine if a point is on one side of the plane or the other

      % Determine if points are inside slab
      dist = sum(xyz .* shape.normal, 1);
      b = abs(dist) <= shape.depth./2;
    end

    function nxyz = normalsXyzInternal(shape, xyz)
      sgn = sign(sum(xyz .* shape.normal, 1));
      nxyz = sgn.*repmat(shape.normal, 1, size(xyz, 2));
    end

    function varargout = intersectAllInternal(shape, vecs)
      % Defer to Strata

      shape = ott.shapes.Strata(shape);
      [varargout{1:nargout}] = shape.intersectAllInternal(vecs);
    end

    function varargout = intersectInternal(shape, vecs)
      % Defer to Strata

      shape = ott.shapes.Strata(shape);
      [varargout{1:nargout}] = shape.intersectInternal(vecs);
    end

    function S = surfInternal(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Converts the shape to a PatchMesh and calls surf.
      %
      % Usage
      %   shape.surfInternal(...)
      %
      % Optional named parameters
      %   - scale (numeric) -- Size of patch.  Default: ``1.0``.
      %
      % Additional named parameters are passed to PatchMesh.surfInternal.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shapes.PatchMesh(shape, 'scale', p.Results.scale);
      S = shape.surfInternal(unmatched{:});
    end

    function shape = scaleInternal(shape, sc)
      shape.depth = shape.depth * sc;
    end
  end

  methods % Getters/setters
    function shape = set.depth(shape, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'depth must be positive numeric scalar');
      shape.depth = val;
    end

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
    function n = get.normal(shape)
      n = shape.rotation(:, 3);
    end

    function bb = get.boundingBox(shape)
      bb = [-Inf, Inf; -Inf; Inf; -shape.depth./2; shape.depth./2];
    end

    function b = get.starShaped(~)
      % Finite depth, so it is star shaped
      b = true;
    end
    function b = get.xySymmetry(~)
      % Origin is at centre of slab
      b = true;
    end
    function q = get.zRotSymmetry(~)
      q = 0;
    end
  end
end
