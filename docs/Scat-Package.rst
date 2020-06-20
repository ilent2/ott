
.. automodule:: +ott.+scat

.. _scat-package:

##############
`scat` Package
##############

The :mod:`+ott.+scat` package contains methods for simulating optical
scattering.
The package is split into multiple sub-packages for the different
scattering methods and a :mod:`+utils` package for mixin and support
classes used for declaring methods.
A summary of the different scattering methods is shown in
:numref:`package-overview-scat`.

.. _package-overview-scat:
.. figure:: images/packageOverview/scat.png
   :alt: Graphical (overview) table of contants for scat package

   Graphical overview of the different scattering methods.
   Different methods deal with different kinds of problems: some
   methods only work well for certain size ranges, while other methods
   require certain types of beams or particles.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   Tmatrix-Package.rst
   Dipole-Package.rst
   Geometric-Package.rst
   Planewave-Package.rst
   Interp-Package.rst
   Shapeforce-Package.rst
   Scat-Utils-Package.rst

