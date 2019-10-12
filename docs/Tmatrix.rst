
#################
`Tmatrix` classes
#################

This page contains information about the ``Tmatrix`` classes currently
implemented in the toolbox.

.. contents::
   :depth: 3
..


Tmatrix
=======

Base class representing T-matrix of a scattering particle or lens. This
class is the base class for all other T-matrix object, you should
inherit from this class when defining your own T-matrix creation
methods. This class doesn't inherit from ``double`` or ``single``,
instead the internal array type can be set at creation allowing the use
of different data types such as ``sparse`` or ``gpuArray``.

TmatrixEbcm
===========

Constructs a T-matrix using extended boundary conditions method.

TmatrixMie
==========

Construct T-matrix from Mie scattering coefficients.

TmatrixPm
=========

Constructs a T-matrix using the point matching method.

TmatrixSmarties
===============

TmatrixDda
==========
