classdef TriangularMesh < ott.shapes.Shape ...
    & ott.shapes.mixin.CoordsCart ...
    & ott.shapes.mixin.VarStarShaped ...
    & ott.shapes.mixin.VarXySymmetry ...
    & ott.shapes.mixin.VarZRotSymmetry ...
    & ott.shapes.mixin.Patch ...
    & ott.shapes.mixin.IntersectMinAll
% Describes a mesh formed by triangular patches.
%
% This class is similar to :class:`PatchMesh` except the patches
% must be triangles (described by three vertices).
%
% Properties
%   - verts     -- 3xN matrix of vertex locations
%   - faces     -- 3xN matrix of vertex indices describing faces
%   - norms     -- 3xN matrix of face normal vectors
%
% Methods
%   - subdivide  -- Add an extra vertex to the centre of each face
%
% Faces vertices should be ordered so normals face outwards for
% volume and inside functions to work correctly.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    verts           % Matrix of vertices in the object
    faces           % Matrix of faces in the object
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    norms              % Matrix of face normal vectors
  end

  methods
    function shape = TriangularMesh(verts, faces, varargin)
      % Construct a new triangular mesh representation
      %
      % TriangularMesh(verts, faces)
      %   verts       3xN matrix of vertex locations
      %   faces       3xN matrix of vertex indices describing faces
      %
      % Faces vertices should be ordered so normals face outwards for
      % volume and inside functions to work correctly.

      shape = shape@ott.shapes.Shape(varargin{:});

      shape.verts = verts;
      shape.faces = faces;
    end

    function shape = subdivide(shape, varargin)
      % Add an extra vertex to the centre of each face
      %
      % Usage
      %   shape = shape.subdivide(...)

      ott.utils.nargoutCheck(shape, nargout);

      v1 = shape.verts(:, shape.faces(1, :));
      v2 = shape.verts(:, shape.faces(2, :));
      v3 = shape.verts(:, shape.faces(3, :));

      v4 = (v1 + v2 + v3)./3;
      idx = (1:size(v4, 2)) + size(shape.verts, 2);

      shape.verts = [shape.verts, v4];

      numFaces = size(shape.faces, 2);
      shape.faces = repmat(shape.faces, 1, 3);

      shape.faces(1, (1:numFaces)+0*numFaces) = idx;
      shape.faces(2, (1:numFaces)+1*numFaces) = idx;
      shape.faces(3, (1:numFaces)+2*numFaces) = idx;
    end

    function [xyz, nxyz, dA] = surfPoints(shape, varargin)
      % Calculate point on surface for surface integration
      %
      % This places one point per face.  It may be better to subdivide
      % faces or use an angular grid depending on the application.
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)

      p = inputParser;
      p.parse(varargin{:});

      % Get raw locations
      v1 = shape.verts(:, shape.faces(1, :));
      v2 = shape.verts(:, shape.faces(2, :));
      v3 = shape.verts(:, shape.faces(3, :));

      xyz = (v1 + v2 + v3)./3;
      nxyz = cross(v1-v2, v3-v2);
      dA = vecnorm(nxyz);
      nxyz = nxyz ./ dA;
    end

    function writeWavefrontObj(shape, filename)
      % Write representation of shape to Wavefront OBJ file
      %
      % Usage
      %   shape.writeWavefrontObj(filename)

      fp = fopen(filename, 'w');

      % Write documentation
      fprintf(fp, '# Shape description generated by OTT\n');

      % Write verts
      for ii = 1:size(shape.verts, 2)
        fprintf(fp, 'v %f %f %f\n', shape.verts(:, ii));
      end

      % Write faces
      for ii = 1:size(shape.faces, 2)
        fprintf(fp, 'f ');
        fprintf(fp, '%d ', shape.faces(:, ii));
        fprintf(fp, '\n');
      end

      % Close file
      fclose(fp);
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz)
      % Determine if Cartesian point is inside the shape

      % Using a third-party function for insidexyz
      b = ott.utils.inpolyhedron(shape.faces.', shape.verts.', xyz.');
    end

    function norms = normalsXyzInternal(shape, xyz)
      % Determine normals for point
      %
      % This is a similar procedure to intersection testing except the
      % vector direction is the normal direction.

      % Calculate normals (3xN)
      N = shape.norms;

      % Reshape points array
      Q1 = reshape(xyz, 3, 1, []);

      % Ensure size of D and N match
      N = repmat(N, 1, 1, size(Q1, 3));
      Q1 = repmat(Q1, 1, size(N, 2), 1);
      
      % Get vertex coordinates (3xN)
      P1 = repmat(shape.verts(:, shape.faces(1, :)), 1, 1, size(Q1, 3));
      P2 = repmat(shape.verts(:, shape.faces(2, :)), 1, 1, size(Q1, 3));
      P3 = repmat(shape.verts(:, shape.faces(3, :)), 1, 1, size(Q1, 3));

      % Calculate intersection points (3xNxM)
      dist = dot(P1 - Q1, N);
      P = Q1 + dist.*N;

      % Determine which points are inside triangles
      % Add a tolerance to make faces overlap slightly
      tol = -1.0e-2;
      I = dot(cross(P2-P1, P-P1), N) >= tol & ...
          dot(cross(P3-P2, P-P2), N) >= tol & ...
          dot(cross(P1-P3, P-P3), N) >= tol;
        
      % Discard points that don't intersect
      dist(:, ~I) = nan;
      N(:, ~I) = nan;

      % Select only nearest face
      [~, idx] = min(abs(dist), [], 2);

      % Select outputs
      norms = N(:, idx);
    end

    function [P, N, dist] = intersectAllInternal(shape, vecs)
      % Find all face intersections, returns a 3xNxM matrix

      % Reshape points (3x1xM)
      Q1 = reshape(vecs.origin, 3, 1, []);
      D = reshape(vecs.direction, 3, 1, []);

      % Calculate normals (3xN)
      N = shape.norms;
      
      % Ensure size of D and N match
      N = repmat(N, 1, 1, size(D, 3));
      D = repmat(D, 1, size(N, 2), 1);
      Q1 = repmat(Q1, 1, size(N, 2), 1);

      % Get vertex coordinates (3xN)
      P1 = repmat(shape.verts(:, shape.faces(1, :)), 1, 1, size(D, 3));
      P2 = repmat(shape.verts(:, shape.faces(2, :)), 1, 1, size(D, 3));
      P3 = repmat(shape.verts(:, shape.faces(3, :)), 1, 1, size(D, 3));

      % Calculate intersection points (3xNxM)
      dist = dot(P1 - Q1, N)./dot(D, N);
      P = Q1 + dist.*D;

      % Determine which points are inside triangles
      I = dot(cross(P2-P1, P-P1), N) >= 0 & ...
          dot(cross(P3-P2, P-P2), N) >= 0 & ...
          dot(cross(P1-P3, P-P3), N) >= 0;

      % Remove intersections in the opposite direction
      found = dist >= 0 & I;
      dist(~found) = nan;
      P(:, ~found) = nan;
      N(:, ~found) = nan;

    end

    function varargout = intersectInternal(shape, varargin)
      % Disambiguate
      [varargout{1:nargout}] = intersectInternal@ ...
          ott.shapes.mixin.IntersectMinAll(shape, varargin{:});
    end
  end

  methods % Getters/setters
    function shape = set.verts(shape, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'verts must be 3xN numeric matrix');
      shape.verts = val;
    end

    function shape = set.faces(shape, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'faces must be 3xN numeric matrix');

      if max(shape.faces(:)) > numel(shape.verts)
        error('faces matrix refers to non-existent vertices');
      end

      shape.faces = val;
    end

    function r = get.maxRadius(shape)
      r = max(vecnorm(shape.verts));
    end

    function totalVolume = get.volume(shape)
      % Calculate the volume of the shape
      %
      % This function is based on a surface triangulation volume code
      % version 1.0.0.0 (1.43 KB) by Krishnan Suresh
      % matlabcentral/fileexchange/26982-volume-of-a-surface-triangulation
      % See tplicenses/stl_KrishnanSuresh.txt for information about licensing.

      p = shape.verts;
      t = shape.faces;

      % Compute the vectors d13 and d12
      d13= [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:))); ...
          (p(3,t(2,:))-p(3,t(3,:)))];
      d12= [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); ...
          (p(3,t(1,:))-p(3,t(2,:)))];

      cr = cross(d13,d12,1);%cross-product (vectorized)

      area = 0.5*sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);% Area of each triangle
      totalArea = sum(area);

      crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
      zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
      nz = -cr(3,:)./crNorm;% z component of normal for each triangle

      volume = area.*zMean.*nz; % contribution of each triangle
      totalVolume = sum(volume);%divergence theorem
    end

    function bb = get.boundingBox(shape)
      bb = [min(shape.verts(1, :)), max(shape.verts(1, :));
            min(shape.verts(2, :)), max(shape.verts(2, :));
            min(shape.verts(3, :)), max(shape.verts(3, :))];
    end

    function nxyz = get.norms(shape)
      V1 = shape.verts(:, shape.faces(2, :)) ...
          - shape.verts(:, shape.faces(1, :));
      V2 = shape.verts(:, shape.faces(3, :)) ...
          - shape.verts(:, shape.faces(1, :));
      nxyz = cross(V1, V2);
      nxyz = nxyz ./ vecnorm(nxyz);
    end
  end
end
