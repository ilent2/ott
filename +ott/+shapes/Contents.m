% ott.shapes Descriptions of geometric shapes
%
% Base classes
%   Shape          - Abstract base class for shapes
%   TriangularMesh - Base for triangular mesh objects (such as file loaders)
%   StarShape      - Abstract base class for star shaped particles
%   AxisymShape    - Abstract base class for axis-symmetric particles
%
% Simple geometric shapes
%   Cube           - Cube shape
%   Cylinder       - Cylinder shape
%   Ellipsoid      - Ellipsoid shape
%   Sphere         - Sphere shape
%   Superellipsoid - Superellipsoid shape
%   BiconcaveDisc  - Biconcave disc
%
% File loaders
%   StlLoader      - Loads a shape from a STL file
%   WavefrontObj   - Loads a shape from a Wavefront OBJ file
%
% Shape builders
%   AxisymLerp     - A axis-symmetric particle with lerping between points
%
% Collections of shapes
%   Union          - Represents union between two shapes
%
% Copyright 2018-2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.
