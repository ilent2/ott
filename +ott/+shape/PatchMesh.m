classdef PatchMesh < ott.shape.Shape ...
    & ott.shape.mixin.CoordsCart ...
    & ott.shape.mixin.VarStarShaped ...
    & ott.shape.mixin.VarXySymmetry ...
    & ott.shape.mixin.VarZRotSymmetry ...
    & ott.shape.mixin.Patch
% A surface resembling Matlab polygon patches.
%
% This surface casts to :class:`TriangularMesh` for most operations.
%
% Properties
%   - verts     -- (3xN numeric) Array of vertices for forming faces
%   - faces     -- (mxN numeric) Vertex indices for polygons

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    verts
    faces
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
  end

  methods (Static)
    function shape = FromSurfMatrix(X, Y, Z, varargin)
      % Construct a Patch from three surf compatible matrices
      %
      % Usage
      %   shape = FromSurfMatrix(X, Y, Z, ...)

      verts = [X(:), Y(:), Z(:)];
      [verts, ~, idx] = unique(verts, 'rows');
      verts = verts.';
      idx = reshape(idx, size(X));

      numFaces = prod(size(X)-1);
      faces = zeros(4, numFaces);
      for ii = 1:size(X, 2)-1
        for jj = 1:size(X, 1)-1
          faceIdx = idx(jj:jj+1, ii:ii+1);
          faces(:, jj+(ii-1)*(size(X, 1)-1)) = faceIdx([1,3,4,2]);
        end
      end

      shape = ott.shape.PatchMesh(verts, faces, varargin{:});
    end
  end

  methods
    function shape = PatchMesh(verts, faces, varargin)
      % Construct a new Patch mesh representation
      %
      % Usage
      %   shape = PatchMesh(verts, faces, ...)

      shape = shape@ott.shape.Shape(varargin{:});

      shape.verts = verts;
      shape.faces = faces;
    end

    function shape = ott.shape.TriangularMesh(shape)
      % Convert the shape to a TriangularMesh
      %
      % Converts faces to triangles and creates a new Triangular mesh.

      % Convert to triangles
      faces = shape.faces(1:3, :);
      for ii = 4:size(shape.faces, 1)
        faces = [faces, [shape.faces([1, ii-1, ii], :)]];
      end

      shape = ott.shape.TriangularMesh(shape.verts, faces, ...
          'position', shape.position, 'rotation', shape.rotation);
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz)
      % Determine if Cartesian point is inside the shape

      shape = ott.shape.TriangularMesh(shape);
      b = shape.insideXyzInternal(xyz);
    end

    function nxyz = normalsXyzInternal(shape, xyz)
      % Determine normals for point

      shape = ott.shape.TriangularMesh(shape);
      nxyz = shape.normalsXyzInternal(xyz);
    end

    function shape = scaleInternal(shape, sc)
      shape.verts = shape.verts .* sc;
    end
  end

  methods % Getters/setters

    function shape = set.verts(shape, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'verts must be 3xN numeric matrix');
      shape.verts = val;
    end

    function shape = set.faces(shape, val)
      assert(isnumeric(val) && ismatrix(val), ...
          'faces must be mxN numeric matrix');

      if max(shape.faces(:)) > numel(shape.verts)
        error('faces matrix refers to non-existent vertices');
      end

      shape.faces = val;
    end

    function bb = get.boundingBox(shape)
      bb = [min(shape.verts(1, :)), max(shape.verts(1, :));
            min(shape.verts(2, :)), max(shape.verts(2, :));
            min(shape.verts(3, :)), max(shape.verts(3, :))];
    end

    function r = get.maxRadius(shape)
      r = max(vecnorm(shape.verts));
    end

    function v = get.volume(shape)
      shape = ott.shape.TriangularMesh(shape);
      v = shape.volume();
    end
  end
end

