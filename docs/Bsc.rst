
.. _bsc-classes:

#############
`Bsc` classes
#############

This section contains information about the beam shape coefficient classes
(Bsc) currently implemented in the toolbox.
These classes can be used to describe optical tweezers beams in a
basis of vector spherical wave functions.
The classes provide functions for translating beams, visualising
beams and overloads for adding beams.
Most of the core functionality is provided in the base class
:class:`+ott.Bsc`.
Classes inheriting from this class typically only need to define
the beam creation code specific to that type of beam.

.. contents:: Contents
   :depth: 3
   :local:
..


Bsc
===

:class:`+ott.Bsc` is the base class for objects representing
beam shape coefficients (BSC) including 
:class:`BscPmGauss` and :class:`BscPmParaxial`.
The class can also be used directly to package a set of existing
BSC for use with other functions in the toolbox, for example

.. code-block:: matlab

   a = [1; 0; 0]; b = 1i.*a;
   basis = 'incoming';
   type = 'incident';
   beam = ott.Bsc(a, b, basis, type);

would create a new beam with `Nmax = 1` (i.e. 3 coefficients for `a` and `b`)
with the incoming vector spherical wave function basis,
representing a incident beam.
For further information about creating custom beams, see
the :ref:`creating-a-custom-beam` example.

.. autoclass:: +ott.Bsc
   :members: Bsc, GetVisualisationData, visualise,
      visualiseFarfield, visualiseFarfieldSlice,
      visualiseFarfieldSphere, translateXyz, translateZ

BscPlane
========

Representation of a plane wave in VSWF coefficients

.. autoclass:: +ott.BscPlane

BscPmGauss
==========

Provides HG, LG and IG beams using point matching method

.. autoclass:: +ott.BscPmGauss
   :members: BscPmGauss

BscPmParaxial
=============

Calculate representation from farfield/paraxial beam

.. autoclass:: +ott.BscPmParaxial

BscPointmatch
=============

Base class for BSC generated using point matching

.. autoclass:: +ott.BscPointMatch

