ott - Optical Tweezers Toolbox
==============================

The optical tweezers toolbox can be used to calculate optical forces
and torques of particles using the T-matrix formalism in a vector
spherical wave basis.
The toolbox includes codes for calculating T-matrices, beams described
by vector spherical wave functions, functions for calculating forces
and torques, and examples.

Installation
------------

To use the toolbox, download or clone the GitHub repository.
The toolbox consists of a directory, ott-toolbox, which contains
a collection of functions for calculating T-matricies, beam coefficients,
force and torques.
To use the functions in your code, the easiest way is to add the
ott-toolbox directory to your matlab path,

```matlab
addpath('<download-path>/ott/ott-toolbox');
```

if you regularly use the toolbox you might want to add the command to
your [startup.m](https://au.mathworks.com/help/matlab/ref/startup.html?searchHighlight=startup.m) file.
You might also want to add the examples to your path

```matlab
addpath('<download-path>/ott/examples');
```

Getting started
---------------

1. To get started using the toolbox you should first take a look at the
users guide for version 1.2 and the optical tweezers computational
toolbox paper (pre-print).  Both are available on our
[website](https://people.smp.uq.edu.au/TimoNieminen/software.html).

2. Copy the examples to your working directory, and play with
    them. Start with example_gaussian.m.

3. It's best to use length units of the wavelength in the trapping
    medium, usually free-space wavelength/1.33.
    
4. The examples calculate the force and torque efficiencies. These are
    the force and torque per photon, in photon units. To convert to
    SI units:
                force_SI = force_Q * n * P/c
               torque_SI = torque_Q * P/w
    where n is the refractive index of the surrounding medium,
          P is the beam power in watts,
          c is the speed of light in free space,
          w is the angular optical frequency, in radians/s.

Upgrading
---------

* Download the latest release of 1.3, add it to your Matlab path with

```
addpath('ott-1.3/ott-toolbox');
```

* Make sure change warnings are on, i.e. add the following to your code:

```
ott_warning('once');
change_warnings('on');
```

* Run your code and take note of the warnings produced
* Download the latest release of 1.4.  Remove 1.3 from your matlab path
  and add 1.4 to your matlab path.

```
rmpath('ott-1.3/ott-toolbox');
addpath('ott-1.4/ott-toolbox');
```

* Turn off change warnings, since these are now the warnings about
  depreciations/changes to 1.5.

```
ott_warning('off');
ott_warning('on');
ott_warning('once');
change_warnings('off');
```

* For each of the warnings that you noted before when you ran your
  code for 1.3, change the appropriate lines.

  Common changes include:

  * `lg_mode_w0` is now depreciated.  If you are using the output of
    `lg_mode_w0` as input to the `bsc_*` functions, you should now pass
    the `beam_angle` into these functions instead of `w0`.

  * `z_rotation_matrix` can be replaced with `rotz(phi_deg)*roty(theta_deg)`
    where `phi_deg` and `theta_deg` are the azimuthal and axial
    rotations in degrees.

  * There is now only one force/torque calculation function, `forcetorque`.

  * `calc_rotation_matrix` can be replaced with `rotation_matrix`.

* Run your code again, checking the result is correct.

Upcoming release
----------------

* Version 1.5 will introduce an object orientated interface for
  beams and T-matrices.  Nmax will be hidden within the objects,
  and automatic choice of Nmax will be done where possible.

* Version 2 will introduces a focus on simulating particles in
  optical traps rather than just focussing on calculating optical
  forces and torques.  The plan is also to introduce geometric
  optics, Rayleigh particles, and arbitrary T-matrix particles
  calculated using DDA.  The toolbox will be more automated and
  include a graphical user interface.

Licence
-------

The package and its components may be used free-of-charge for research,
teaching, or personal use. If results obtained using the package are
published, the package should be appropriately referenced.
For full details see LICENSE.md.
For use outside the conditions of the license, please contact us.

The package can be refereced by citing the paper describing version
1 of the toolbox

> T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, G. Knöner, A. M. Branczyk, N. R. Heckenberg, and H. Rubinsztein-Dunlop,
> "Optical tweezers computational toolbox",
> [Journal of Optics A 9, S196-S203 (2007)](http://iopscience.iop.org/1464-4258/9/8/S12/)

or by directly citing the toolbox

> T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, I. C. D. Lenton,
> Y. Hu, G. Knoener, A. M. Branczyk,
> "Optical tweezers toolbox", https://github.com/ilent2/ott

Contact us
----------

The best person to contact for inquiries about the toolbox or lincensing
is [Timo Nieminen](mailto:timo@physics.uq.edu.au)

Further Reading
---------------

Papers describing the toolbox

* T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, G. Knoener,
A. M. Branczyk, N. R. Heckenberg, H. Rubinsztein-Dunlop,
"Optical tweezers computational toolbox",
Journal of Optics A 9, S196-S203 (2007)

* T. A. Nieminen, V. L. Y. Loke, G. Knoener, A. M. Branczyk,
"Toolbox for calculation of optical forces and torques",
PIERS Online 3(3), 338-342 (2007)


More about computational modelling of optical tweezers:

* T. A. Nieminen, N. R. Heckenberg, H. Rubinsztein-Dunlop,
"Computational modelling of optical tweezers",
Proc. SPIE 5514, 514-523 (2004)


More about our beam multipole expansion algorithm:

* T. A. Nieminen, H. Rubinsztein-Dunlop, N. R. Heckenberg,
"Multipole expansion of strongly focussed laser beams",
Journal of Quantitative Spectroscopy and Radiative Transfer 79-80,
1005-1017 (2003)

More about our T-matrix algorithm:

* T. A. Nieminen, H. Rubinsztein-Dunlop, N. R. Heckenberg,
"Calculation of the T-matrix: general considerations and
application of the point-matching method",
Journal of Quantitative Spectroscopy and Radiative Transfer 79-80,
1019-1029 (2003)

The multipole rotation matrix algorithm we used:

* C. H. Choi, J. Ivanic, M. S. Gordon, K. Ruedenberg,
"Rapid and stable determination of rotation matrices between
spherical harmonics by direct recursion"
Journal of Chemical Physics 111, 8825-8831 (1999)


The multipole translation algorithm we used:

* G. Videen,
"Light scattering from a sphere near a plane interface",
pp 81-96 in:
F. Moreno and F. Gonzalez (eds),
Light Scattering from Microstructures, LNP 534,
Springer-Verlag, Berlin, 2000


More on optical trapping landscapes:

* A. B. Stilgoe, T. A. Nieminen, G. Knoener, N. R. Heckenberg, H. 
Rubinsztein-Dunlop, "The effect of Mie resonances on trapping in 
optical tweezers", Optics Express, 15039-15051 (2008)

Multi-layer sphere algorithm:

* W. Yang, "Improved recursive algorithm for light scattering by a
 multilayered sphere", Applied Optics 42(9), (2003)

