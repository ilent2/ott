
.. _getting-started:

###############
Getting Started
###############

This section will guide you through installing the toolbox and running
your first toolbox commands either through the graphical user
interface (GUI) or on the command line.

.. contents:: Contents
   :depth: 3
   :local:
..

Installation
============

There are several methods for installing the toolbox.
If you are only interested in using the latest stable release of the toolbox,
the easiest method is via the Matlab Addons Explorer.
Alternatively, you can download a specific release from the
`OTT GitHub page <https://github.com/ilent2/ott>`__.
The following sections will guide you through installing the toolbox
and verifying that it is on the Matlab path.

Installing via Matlab Addons Explorer
-------------------------------------

The easiest method to install the toolbox is using the Matlab Addons
explorer.
Simply launch Matlab and navigate to Home > Addons > Get-Addons and
search for "Optical Tweezers Toolbox".
Then, simply click the "Add from GitHub" button to automatically download
the package and add it to the path. You may need to log-in or create a
Mathworks account to complete this step.

Using a .zip or cloning the repository
--------------------------------------

The latest version of OTT can be downloaded from the
`OTT GitHub page <https://github.com/ilent2/ott>`__.
Either select the "Code" button near the top of the screen and select
your preferred download method, or navigate to the
`release page <https://github.com/ilent2/ott/releases>`__
and select the .zip file for the desired release.
Alternatively, you can clone the git repository directly.
There are a range of online tutorials for getting started with git
and GitHub, for example
https://product.hubspot.com/blog/git-and-github-tutorial-for-beginners.
If you downloaded a .zip file, you will need to extract the files
to somewhere on your computer before proceeding.

Next, you need to tell Matlab where to find OTT.
To do this, simply run

.. code:: matlab

    addpath('/path/to/toolbox/ott');

Replace the path with the path you placed the downloaded toolbox in. The
folder must contain the ``+ott`` directory and the ``docs`` directory.
If you downloaded the latest toolbox from GitHub, the final part of the
pathname will either be the repository path (if you used ``git clone``)
or something like ``ott-master`` (if you downloaded a .zip file). The above
line can be added to the start of each of your files or for a more
permanent solution you can add it to the `Matlab startup
script <https://au.mathworks.com/help/matlab/ref/startup.html>`__.

Verifying the installation
--------------------------

To verify that the toolbox can be found by Matlab, type

.. code:: matlab

  help ott

into the Matlab command window, which should display the contents of
``+ott/Contents.m`` if everything installed correctly.
If you have multiple versions of the
toolbox installed and you wish to check which version is currently
active, you can type

.. code:: matlab

  what ott

Using the graphical user interface
==================================

The toolbox includes a graphical user interface for many of the most
common tasks.  The user interface applications are located in the
``ott.ui`` sub-package.  The easiest way to launch the graphical user
interfaces is via the Launcher app.  To start the Launcher, simply
run the following on the Matlab command line (after the toolbox has
been installed):

.. code:: matlab

  ott.ui.Launcher

.. todo:: Update the how to cite text in the Launcher (and finish launcher)

Using the toolbox functions
===========================

The toolbox contains a collection of functions and classes for performing
a range of optical tweezers related simulation tasks.
Different tasks, such as simulating drag or describing geometric shapes,
are split into different sub-packages.
Details about these sub-packages can be found in the :ref:`reference`
section or a short list will be printed to the screen when you run
``help ott``.
Within each package are either functions, classes, or other sub-packages
which further group functions/classes based on their purpose.
For example, the ``beam`` sub-package contains several classes for
generating beams.
In order to create a new beam instance you can either use the class
constructor or a suitable static method (if provided).
For example, to create a new Gaussian beam you would call either the
``Gaussian`` or ``FromNa`` method of the Gaussian class, for example

.. code:: matlab

  beam1 = ott.beam.Gaussian()
  beam2 = ott.beam.Gaussian.FromNa(0.9)

In both cases you need to prefix the class name with the package name.
If you intend to use a range of methods from one package, it is possible
to import that package using

.. code:: matlab

  import ott.beam.*

Except for static functions (such as the ``FromNa`` method above),
most class functions cannot be called without an instance of the class.
For example, all ``ott.beam.Beam`` object implement a function called
``efield`` which calculates the electric field around the coordinate
origin.
In order to use this function you need to first construct a valid
beam object (for example, using the ``Gaussian`` or ``Gaussian.FromNa``
methods above).
The following example shows how to create a new Gaussian beam and
calculate the field near the origin with the ``efield`` method.

.. code:: matlab

  beam = ott.beam.Gaussian()
  E = beam.efield([0;0;1e-6])    % Calculate field near origin

The ``examples`` directory includes multiple examples demonstrating
various features of the toolbox.
To get started writing your own code, we recommend that you start by
working through the examples and reading the :ref:`examples` section of
this manual.
To get help on a specific method or class, you can either type
``help <name of method>`` or lookup the method/class in the
:ref:`reference` section of this manual.
For further information on using Matlab packages and classes, refer to
the `Mathworks OOP documentation <https://au.mathworks.com/help/matlab/object-oriented-programming.html>`__.

