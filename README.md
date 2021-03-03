ott - Optical Tweezers Toolbox
==============================

[![Documentation Status](https://readthedocs.org/projects/ott/badge/?version=latest)](https://ott.readthedocs.io/en/latest/?badge=latest)
[![View Optical Tweezers Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://au.mathworks.com/matlabcentral/fileexchange/73541-optical-tweezers-toolbox)

The Optical Tweezers Toolbox is a collection of methods for modelling
optically trapped particles.
The toolbox includes methods for modelling
tightly focussed laser beams, optical scattering, fluidic drag, and
dynamics of optically trapped particles.
The toolbox can be used to calculate the optical forces and torques
acting on spherical and non-spherical particles in a laser beam.
This version of the toolbox includes methods for calculating fluidic
drag and Brownian forces, as well as examples for basic dynamics
simulations.

Installation and usage
----------------------

There are several methods for installing the toolbox.
For the latest stable release, the easiest installation method
is to launch Matlab and navigate to Home -> Addons -> Get-Addons
and search for "Optical Tweezers Toolbox".  Then, simply click the
"Add from GitHub" button to automatically download the package and
add it to the path.
Alternatively, you can download the toolbox directly from the
[GitHub repository](https://github.com/ilent2/ott) or select a
specific release; if using this method you will need to add the
toolbox directory to Matlab's path using

```matlab
addpath('<download-path>/ott');
help ott   % Test that ott was found, should display ott Contents page
```

Regardless of the method you acquired OTT with: once installed you
should be able to access the toolbox functions and classes contained
in the `ott.*` package, for example

```matlab
beam = ott.beam.Gaussian();
```

or for the graphical user interface

```matlab
ott.ui.Launcher
```

More detailed instructions can be found in the Getting Started
section of the [documentation](https://ott.readthedocs.io/).

Dependencies
------------

The toolbox should be compatible with **Matlab 2018a or newer**.
Previous versions of the toolbox were compatible with
[GNU Octave](https://www.gnu.org/software/octave/), however the current
version does not appear to be compatible in the most recent test.

Licence
-------

Except where otherwise noted, this toolbox is made available under the
Creative Commons Attribution-NonCommercial 4.0 License.
For full details see LICENSE.md.
For use outside the conditions of the license, please contact us.
The toolbox includes some third-party components, information about
these components can be found in the documentation and corresponding
file in the thirdparty directory.

If you use the toolbox in academic work, please cite it as follows

TODO

or using the following bibtex entry

TODO

Contact us
----------

The best person to contact for inquiries about the toolbox or lincensing
is [Isaac Lenton](mailto:uqilento@uq.edu.au).

File listing
------------

```
README.md     - Overview of the toolbox (this file)
LICENSE.md    - License information
AUTHORS.md    - List of contributors and info about copyright
CHANGES.md    - Overview of changes to the toolbox
thirdparty/   - Third party licenses (multiple files)
examples/     - Example files showing different toolbox features
tests/        - Unit tests to verify toolbox features function correctly
docs/         - Sphinx documentation for the project
+ott/         - The toolbox (formatted as a Matlab pacckage)
```

The +ott package, as well as tests/, examples/, docs/ and thirdparty/
directories and sub-directories contain Contents.m files which list
the files and packages in each directory.
These files can be viewed in Matlab by typing `help ott`
or `help ott.subpackage`.
