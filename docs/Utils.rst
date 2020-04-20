
.. automodule:: +ott.+utils

.. _utils-package:

###############
`utils` Package
###############

This package contains functions used by other methods in the toolbox.
Most of these functions are from version 1 of the toolbox with some
minor modifications and bug fixes.

.. warning:: These functions are likely to move in future releases and are
   not very well documented.

.. contents::
   :depth: 1
   :local:
..

Special functions
=================

.. autofunction:: sbesselh1
.. autofunction:: sbesselh2
.. autofunction:: sbesselh
.. autofunction:: sbesselj

.. autofunction:: spharm
.. autofunction:: vsh
.. autofunction:: vswf
.. autofunction:: vswfcart

Coordinate transformations
==========================

.. autofunction:: rtp2xyz
.. autofunction:: rtpv2xyzv

.. autofunction:: xyz2rtp
.. autofunction:: xyzv2rtpv

Translations and rotations
==========================

.. autofunction:: translate_z
.. autofunction:: rotx
.. autofunction:: roty
.. autofunction:: rotz

.. autofunction:: rotation_matrix
.. autofunction:: wigner_rotation_matrix

Helper functions
================

.. autofunction:: matchsize
.. autofunction:: threewide
.. autofunction:: iseven
.. autofunction:: isodd
.. autofunction:: rotate_3x3tensor

.. autofunction:: ka2nmax
.. autofunction:: nmax2ka

=================

Geometry functions
==================

.. autofunction:: angulargrid
.. autofunction:: perpcomponent
.. autofunction:: inpolyhedron

Unclassified
============

.. todo:: These functions should be moved to other categories

.. autofunction:: paraxial_transformation_matrix
.. autofunction:: paraxial_beam_waist
.. autofunction:: lgmode
.. autofunction:: legendrerow
.. autofunction:: laguerre
.. autofunction:: incecoefficients
.. autofunction:: hgmode
.. autofunction:: hermite
.. autofunction:: emField
.. autofunction:: combined_index

Polarizability calculation
==========================

.. automodule:: +ott.+utils.+polarizability

.. autofunction:: FCD
.. autofunction:: LDR
.. autofunction:: CM

