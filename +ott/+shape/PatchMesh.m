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

  properties (SetAccess=protected)
    verts         % (3xN numeric) Array of vertices for forming faces
    faces         % (mxN numeric) Vertex indices for polygons
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
      %
      % Parameters
      %   - verts (3xN numeric) -- Array of vertices for forming faces
      %   - faces (mxN numeric) -- Vertex indices for polygons

      shape = shape@ott.shape.Shape(varargin{:});
      shape = shape.setData(verts, faces);
    end
    
    function shape = setData(shape, verts, faces)
      % Set the vertex and face data for the shape
      %
      % Usage
      %   shape = shape.setData(verts, faces)
      %
      % Parameters
      %   - verts (3xN numeric) -- Array of vertices for forming faces
      %   - faces (mxO numeric) -- Vertex indices for polygons
      
      assert(isnumeric(faces) && ismatrix(faces), ...
          'faces must be mxO numeric matrix');
      assert(isnumeric(verts) && ismatrix(verts) && size(verts, 1) == 3, ...
          'verts must be a 3xN numeric matrix');
      assert(min(faces(:)) >= 1 && max(faces(:)) <= size(verts, 2), ...
        	'faces matrix refers to non-existent vertices');
        
      shape.faces = faces;
      shape.verts = verts;
      
    end

    function shape = ott.shape.TriangularMesh(shape)
      % Convert the shape to a TriangularMesh
      %
      % Converts faces to triangles and creates a new Triangular mesh.

      numFaces = size(shape.faces, 2);
      trisPerFace = size(shape.faces, 1)-2;
      
      % Convert to triangles
      ofaces = zeros(3, numFaces*trisPerFace);
      for ii = 1:trisPerFace
        ofaces(:, (1:numFaces) + (ii-1)*numFaces) = ...
            shape.faces([1, ii+1, ii+2], :);
      end

      shape = ott.shape.TriangularMesh(shape.verts, ofaces, ...
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

