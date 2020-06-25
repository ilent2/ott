
.. automodule:: +ott.+beam.+medium

.. _medium-package:

################
`medium` Package
################

This package declares classes for defining optical mediums and materials.
The toolbox uses four different types of classes to declare optical
media:

   * :ref:`mediums` combine a vacuum, material and an optical frequency and
     are used by the beam objects.

   * :ref:`vacuums` declare the relative units to use for beam properties,
     units must have SI dimensions (but can be scaled to have different
     units).

   * :ref:`materials` declare dimensionless material properties as ratios
     of the vacuum medium (i.e. materials are declared with relative
     properties).

   * :ref:`relative-material` defines the ratio between mediums/materials
     and is mainly used by the scattering methods.

To create a beam you will need to specify a medium, a material or a vacuum.
If a material or vacuum is specified, an appropriate medium is constructed
using the defaults for frequency/vacuum.
Mediums can be created without explicitly setting a optical frequency or
vacuum, in this case the default optical frequency and vacuum are used.
These defaults can be changed by modifying the static properties in
the :class:`Medium` class (see :meth:`Medium.DefaultVacuum` and
:meth:`Medium.DefaultFrequency`).

To create a scattering method, you will typically specify a relative
material.  A relative material can be created from two
material/mediums/vacuums using::

   relMaterial = material1 ./ material2;

A range of default materials are declared as static properties of the
:class:`Generic` class, for example::

   water = ott.beam.medium.Generic.Water;

These properties may be useful defaults for your simulations
but you should check that these values are appropriate for your
specific use case.

.. contents:: Contents
   :depth: 3
   :local:

.. _mediums:

Medium
======

.. autoclass:: Medium
   :members: Medium, FromWavelength, DefaultVacuum, DefaultFrequency

.. _vacuums:

Vacuum
======

.. autoclass:: Vacuum
   :members: Vacuum

.. _materials:

Materials
=========

Material (abstract)
-------------------

.. autoclass:: Material
   :members: Material

Dielectric
----------

.. autoclass:: Dielectric
   :members: Dielectric, FromIndex

Arbitrary
---------

.. autoclass:: Arbitrary
   :members: Arbitrary

Generic
-------

.. autoclass:: Generic

.. _relative-material:

Relative material
-----------------

.. autoclass:: Relative
   :members: Relative

