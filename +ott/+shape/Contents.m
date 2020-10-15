% ott.shape Descriptions of geometric shapes
%
% This package provides a collection of simple geometric shapes, methods
% for building arbitrary geometric shapes from functions/points, and
% methods to load geometry from files.
%
% In addition to the simple geometric shapes, some of the shape builder
% classes define static methods for building other commonly used shapes.
% See examples/packageOverview/shapes.m for an overview of available shapes.
%
% Most scattering methods involve surface integrals, volume integrals
% or calculation of surface normals; these classes describe shapes with
% these specific quantities in mind.  The complexity of various
% scattering simulations can often be reduced when the particle is
% star shaped, mirror symmetric or rotational symmetric.
% Consequently, all shapes have methods for querying these properties.
%
% Base class
%   Shape             - Base class for all OTT shapes
%
% Simple geometric shapes
%   Cube              - Cube shape
%   RectangularPrism  - Rectangular based prism
%   Cylinder          - Cylinder shape
%   Ellipsoid         - Ellipsoid shape
%   Sphere            - Sphere shape
%   Superellipsoid    - Superellipsoid shape
%   Plane             - Plane with infinite extent
%   Strata            - Collection of parallel planes
%   Slab              - Parallel planes forming a slab
%   Empty             - Shape with no volume
%
% Shape builders
%   PatchMesh      - Mesh described by vertices and polygon faces
%   TriangularMesh - Mesh described by vertices and triangular faces
%   AxisymInterp   - Axis-symmetric particle with discrete points
%   AxisymFunc     - Axis-symmetric particle with parametric function
%
% File loaders
%   StlLoader      - Loads a shape from a STL file
%   ObjLoader      - Loads a shape from a Wavefront OBJ file
%
% Collections of shapes
%   Set            - (Abstract) Collection of shapes
%   Union          - Represents union between two shapes (| operator)
%   Intersection   - Intersection of multiple shapes (& operator)
%   Inverse        - Shape where inside/outside are flipped (~ operator)
%
% Sub-packages
%   +mixin         - Contains helpers classes for declaring shapes
%
% Copyright 2018-2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.
