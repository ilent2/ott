
################
`shapes` Package
################

This page provides an overview of the shapes currently in the toolbox.

.. contents::
   :depth: 3
..


Base classes
============

Shape - Shape abstract class for optical tweezers toolbox shapes
----------------------------------------------------------------

AxisymShape - AxisymShape abstract class for axisymetric particles
------------------------------------------------------------------

StarShape - StarShape abstract class for star shaped particles
--------------------------------------------------------------

Geometric shapes
================

Cube
----

A simple cube shape

Cylinder
--------

A simple cylinder shape

Ellipsoid
---------

A simple ellipsoid shape

Sphere
------

A simple sphere shape

Superellipsoid
--------------

A simple superellipsoid shape

Procedural shapes
=================

AxisymLerp
----------

A axisymmetric particle with lerping between points.

TriangularMesh
--------------

Base class for triangular mesh objects (such as file loaders). Can also
be called directly with a list of vertices and faces.

File loaders
============

StlLoader
---------

Load a shape from a STL file.

WavefrontObj
------------

Load a shape from a Wavefront OBJ file.
