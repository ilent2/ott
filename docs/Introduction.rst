
.. _introduction:

############
Introduction
############

Welcome to the documentation for the
`Optical Tweezers Toolbox (OTT) <https://github.com/ilent2/ott>`__.
This section provides a brief overview of the toolbox.
Subsequent sections will guide you through installing the toolbox, and
running the basic examples.
The :ref:`reference` section is automatically generated from
the toolbox source code and provides a list of all the major toolbox
functions and usage.
The final section provide additional reference material for toolbox
concepts.

.. contents:: Contents
   :depth: 3
   :local:
..

About
=====

The Optical Tweezers Toolbox is a collection of methods for modelling
optically trapped particles.
The toolbox includes methods for modelling
tightly focussed laser beams, optical scattering, fluidic drag, and
dynamics of optically trapped particles.
The toolbox can be used to calculate the optical forces and torques
acting on spherical and non-spherical particles in a laser beam.
In addition, the toolbox includes methods for calculating fluidic
drag and Brownian forces, as well as examples for basic dynamics
simulations.

This is the second major release of the optical tweezers toolbox.
Since version 1, the toolbox has been substantially re-written in order
to focus on describing *beams* and *particles* instead of specific
numerical methods.
This version of the toolbox implements a more modular design using
modules and objects.
The new modular design is intended to make it easier to integrate other
optical force calculation methods such as geometric optics and
the discrete dipole approximation in a upcoming release.
The other major change from version 1 is the inclusion of methods for
calculating non-optical forces and simulating the dynamics of particles.

Previous releases can be downloaded from the
`OTT GitHub page <https://github.com/ilent2/ott>`__.

Dependencies
============

The toolbox should be compatible with **Matlab 2018a or newer**.
Previous versions of the toolbox were compatible with
`GNU Octave <https://www.gnu.org/software/octave/>`__, however the current
version does not appear to be compatible in our most recent test.

Licence
=======

Except where otherwise noted, this toolbox is made available under the
Creative Commons Attribution-NonCommercial 4.0 License.
For full details see the file ``LICENSE.md``.
For use outside the conditions of the license, please contact us.
The toolbox includes some third-party components, information about
these components can be found in the documentation and corresponding
file in the ``thirdparty`` directory.

If you use the toolbox in academic work, please cite it as follows

.. todo:: Add how to cite, and in README too

or using the following bibtex entry

.. todo:: Add how to cite, and in README too

Contributing
============

If you would like to contribute a feature, report a bug or request we
add something to the toolbox, the easiest way is by `creating a new
issue on the OTT GitHub
page <https://github.com/ilent2/ott/issues>`__.

If you have code you would like to submit, fork the repository, add the
code and open a new issue. This method is preferable to pasting the code
in the issue or sending it to us via email since your contribution
details will remain attached to the commit you send (tracking
authorship).

Contact us
==========

The best person to contact for inquiries about the toolbox or licensing
is `Isaac Lenton <mailto:uqilento@uq.edu.au>`__.

