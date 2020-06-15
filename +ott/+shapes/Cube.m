classdef Cube < ott.shapes.Shape ...
    & ott.shapes.mixin.CoordsCart ...
    & ott.shapes.mixin.Patch
% Simple geometric cube.
% Inherits from :class:`Shape`.
%
% Properties
%   - width        -- Width of the cube
%
% Additional properties inherited from base.

% Copyright 2018-2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    width         % Width of cube
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    starShaped         % True if the particle is star-shaped
    xySymmetry         % True if the particle is xy-plane mirror symmetric
    zRotSymmetry       % z-axis rotational symmetry of particle
    verts              % Vertices describing cube
    faces              % Faces describing cube
  end

  methods
    function shape = Cube(varargin)
      % Construct a cube.
      %
      % Usage
      %   shape = Cube(width, ...)
      %   Parameters can be passed as named arguments.
      %
      % Additional parameters are passed to base.

      p = inputParser;
      p.addOptional('width', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shapes.Shape(unmatched{:});
      shape.width = p.Results.width;
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz)
      assert(size(xyz, 1) == 3, 'xyz must be 3xN matrix');
      b = all(abs(xyz) <= shape.width./2, 1);
    end

    function n = normalsXyzInternal(shape, xyz)

      assert(size(xyz, 1) == 3, 'xyz must be 3xN matrix');

      % Determine which dimensions are inside
      tol = 1.0e-3.*shape.width;
      insideDims = abs(xyz) < (shape.width./2 + tol);

      % Calculate normals
      n = sign(xyz) .* double(~inside);

      % Normalise vectors (creates nans for interior points)
      n = n ./ vecnorm(n);

    end
  end

  methods % Getters/setters
    function shape = set.width(shape, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0, ...
          'width must be positive numeric scalar');
      shape.width = val;
    end

    function r = get.maxRadius(shape)
      r = sqrt(3*(shape.width/2).^2);
    end
    function shape = set.maxRadius(shape, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0, ...
          'maxRadius must be positive numeric scalar');
      shape.width = sqrt(r.^2./3).*2;
    end

    function v = get.volume(shape)
      v = shape.width.^3;
    end
    function shape = set.volume(shape, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0, ...
          'volume must be positive numeric scalar');
      shape.width = val.^(1/3);
    end

    function bb = get.boundingBox(shape)
      bb = [-1, 1; -1, 1; -1, 1].*shape.width./2;
    end

    function b = get.starShaped(shape)
      b = true;
    end
    function b = get.xySymmetry(shape)
      b = true;
    end
    function q = get.zRotSymmetry(shape)
      q = 4;
    end

    function verts = get.verts(shape)
      X = [1, 1, -1, -1, 1, 1, -1, -1];
      Y = [-1, 1, 1, -1, -1, 1, 1, -1];
      Z = [-1, -1, -1, -1, 1, 1, 1, 1];

      verts = [X(:), Y(:), Z(:)].';
    end

    function faces = get.faces(shape)
      faces = [1, 2, 6, 5; 2, 3, 7, 6; 3, 4, 8, 7; 4, 1, 5, 8; ...
          2, 1, 4, 3; 5, 6, 7, 8].';
    end
  end
end
