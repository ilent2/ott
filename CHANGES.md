Change history
==============

Version 1.5
-----------

* Makes the toolbox a package
* Creates a utils sub-package for non-essential functions, these functions may
    be more susceptable to change in future releases.
* Introduces Tmatrix and Bsc classes to manage/generate beams and T-matrices
* Introduces a simple method of calculating T-matrices (Tmatrix.simple())

Version 1.4
-----------

* Removed bloat from the toolbox, some functions are now gone
* Removed the depreciation/change warnings added in 1.3.1
* Replaced bsc_bessel_farfield.m with bsc_bessel.m
* Added a new example for simulating SLM patterns
* `bsc_*` functions now accept a beam angle instead of a beam waist

Version 1.3.1
-------------

* Added function name change warnings and depreciation warnings.
* Added functions from latest in-house development version of the toolbox.
* Various improvements and bug fixes
* Added `ott_warning` and `change_warnings` functions to make warnings
    less obtrusive.  Internal calls to toolbox functions now suppress
    many warning messages.

Version 1.3
-----------

Bugs: Fixed various bugs related to the beam code and vswf. These bugs were minor and 
would only have effects on a small range of calculations.

Added: Axisymmetric extended boundary condition codes. One of the codes can calculate 
and arbitrary perimeter with associated normals and surface area elements and the other
calculates the EBCM T-matrix for that particular perimeter. axisym_boundarypoints.m and
tmatrix_ebcm_axisym.m. By default tmatrix_ebcm_axisym.m calls axisym_boundarypoints.m 
for the desired boundary.

Version 1.2
-----------

Added: multilayered sphere code tmatrix_mie_layered.m. Definitely works for two layers.
Change: added vector m functionality in spharm.m, vsh.m and vswf.m. These are still 
compatible with previous version of the toolbox

Change: added proper LG mode support in
bsc_pointmatch_farfield.m. Previously, only LGpl modes
with p=0 could be used.

Change: re-wrote the example_x.m files to be more computationally
efficient. example_landscape.m and example_cube.m now contain
executable code producing results comparable to Nieminen et al.
2007.

Change: made many functions far more vector friendly. Files made more vectory are:
spharm.m -- in addition the matlab legendre function was replaced by my own one.
sbessel(x).m
vsh.m
vswf.m

Change: changed the readme.txt and readme.m to reflect the changes made to the toolbox
and the addition of a user guide.


Version 1.1
-----------

Bug: different power normalisations used in force_z.m,
forcetorque.m, and elesewhere. force_z.m and forcetorque.m
fixed to conform to the standard normalisation
P = sum( abs(a).^2 + abs(b).^2 ). Optical forces and force
efficiencies calculated with the wrong normalisation are
4 times too low.

Bug: example_gaussian.m and example_lg.m used z(nz) instead
of zeq to determine Nmax before doing radial force calculation.
Results numerically correct, but code runs slower. Now fixed.

Change: example_gaussian.m now uses numerical aperture to
determine w0.

Change: example_lg.m now uses beam convergence angle and beam
mode to determine w0.

Version 1.0
-----------
Initial release.

