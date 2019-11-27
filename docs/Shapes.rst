
.. automodule:: +ott.+shapes

.. _shapes-package:

################
`shapes` Package
################

This section provides an overview of the shapes currently in the toolbox.

.. contents::
   :depth: 1
   :local:
..


Base classes
============

.. autoclass:: Shape
.. autoclass:: AxisymShape
.. autoclass:: StarShape

Geometric shapes
================

.. autoclass:: Cube
.. autoclass:: RectangularPrism
.. autoclass:: Cylinder
.. autoclass:: Ellipsoid
.. autoclass:: Sphere
.. autoclass:: Superellipsoid

Sets of shapes
==============

These classes can be used to create shapes by combining simple geometric
shapes or other shape objects.
For instance, the union class can be used to create a union of two
spheres::

   shape1 = ott.shapes.Sphere(1.0, [0, 0, -2]);
   shape2 = ott.shapes.Sphere(1.0, [0, 0,  2]);
   union = ott.shapes.Union([shape1, shape2]);

.. autoclass:: Union
   :members: Union

.. todo:: We will probably add other sets in future including
   differences or exclusions.

Procedural shapes
=================

.. autoclass:: AxisymLerp
.. autoclass:: TriangularMesh

File loaders
============

.. autoclass:: StlLoader
.. autoclass:: WavefrontObj

