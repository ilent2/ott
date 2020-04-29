% ott.shapes Descriptions of geometric shapes
%
% Abstract base classes
%   Shape          - Shape abstract class for optical tweezers toolbox shapes
%   ShapeCart      - Shape abstract class with helpers for Cartesian coords.
%   ShapeSph       - Shape abstract class with helpers for Spherical coords.
%   TriangularMesh - Base for triangular mesh objects (such as file loaders)
%   StarShape      - Abstract base class for star shaped particles
%   AxisymShape    - Abstract base class for axis-symmetric particles
%   Extrusion      - Abstract class for shapes described by a height function
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
