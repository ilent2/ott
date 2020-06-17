classdef Patch < ott.shapes.mixin.IntersectTriMesh
% Shapes describes by lists of vertices and faces.
% Inherits from :class:`IntersectTriMesh`.
%
% Abstract properties
%   - verts     -- 3xN list of N vertices
%   - faces     -- mxN list of patches
%
% Methods
%   - surfInternal -- Generate a surface using the patch function
%   - surfPoints   -- Cast to TriangularMesh and call surfPoints
%   - intersectInternal -- Casts to TriangularMesh
%   - intersectAllInternal -- Casts to TriangularMesh
%
% Supported casts
%   - TriangularMesh    -- Inherited (uses PatchMesh)
%   - PatchMesh

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Abstract)
    verts
    faces
    position
    rotation
  end

  methods
    function [xyz, nxyz, dA] = surfPoints(shape, varargin)
      % Cast to TriangularMesh and call surfPoints
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)

      shape = ott.shapes.TriangularMesh(shape);
      [xyz, nxyz, dA] = shape.surfPoints(varargin{:});
    end

    function shape = ott.shapes.PatchMesh(shape)
      % Convert the shape to a PatchMesh

      shape = ott.shapes.PatchMesh(shape.verts, shape.faces, ...
          'position', shape.position, 'rotation', shape.rotation);
    end
  end

  methods (Hidden)
    function S = surfInternal(shape)
      % Construct surface struture
      S = struct('Vertices', shape.verts.', ...
          'Faces', shape.faces.', 'FaceColor', 'flat', ...
          'FaceVertexCData', [0.9290 0.6940 0.1250]);
    end
  end
end
