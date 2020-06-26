% ott.shapes.mixin Mixin classes for use with shapes.
%
% Coordinates
%   CoordsCart    - Cartesian shape description helper
%   CoordsSph     - Spherical shape description helper
%
% Properties
%   VarStarShaped     - Variable star shaped property
%   VarZRotSymmetry   - Variable z rotational symmetry property
%   VarXySymmetry     - Variable xy mirror symmetry property
%   NumericalVolume   - Calculate volume numerically
%   InfVolume         - Shapes with infinite volume
%
% Abstract properties
%   IsSphereAbsProp   - Adds a isSphere abstract property
%
% Surfaces
%   IsosurfSurf       - Use isosurf for surf
%   IsosurfSurfPoints - use isosurf for surfPoints
%   NoSurfPoints      - Shapes where surfPoints doesn't make sense
%
% Intersection
%   IntersectMinAll   - Declare intersect method that uses intersectAll
%   IntersectRayMarch - Intersection using ray marching
%   IntersectTriMesh  - Intersection by casting to a TriangularMesh
%
% Shape types
%   Patch          - Shape described by a discrete set of vertices
%   StarShape      - Star shaped particles
%   AxisymShape    - Axis-symmetric particles
%   AxisymStarShape  - Combines AxisymShape with StarShape
%
% Copyright 2018-2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

