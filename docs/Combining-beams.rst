
###############
Combining Beams
###############

The toolbox can calculate the force on a particle from multiple beams
either coherently or incoherently. This page describes the different
beam combination methods used in the ``multiple_beams.m`` example. For
both coherent and incoherent beams, the order which beams are combined
and translated can affect accuracy and computation time. This page
assumes the beams can be exactly represented, as is the case for focused
beams which pass through a finite aperture. For plane waves and other
beams requiring an infinite VSWF expansion, it is better to calculate
the beams at the new location.

For this example, we calculate the forces on a spherical partice, the
T-matrix for this particle is given by

.. code:: matlab

    % Wavelength in medium/vacuum [m]
    wavelength = 1064.0e-9;

    T = ott.Tmatrix.simple('sphere', wavelength, 'index_medium', 1.0, ...
        'index_particle', 1.2, 'wavelength0', wavelength);

We combine two copies of the same beam displaced along the x axis:

.. code:: matlab

    beam = ott.BscPmGauss('polarisation', [1 1i], 'angle_deg', 50, ...
        'index_medium', 1.0, 'wavelength0', wavelength, 'power', 0.5);

    % Displacement of beams [wavelength_medium]
    displacement = 0.2*wavelength;

The original beam is centered around the coordinate origin. When we use
these beams we need to translate them to the particle origin. For
example, if the particle is displaced a distance ``x`` from the centre
of the two beams, we can translate the beams to this location using by
translating each beam separately:

.. code:: matlab

    beam1 = beam.translateXyz([x+displacement; 0; 0]);
    beam2 = beam.translateXyz([x-displacement; 0; 0]);

Combining coherent beams
========================

For coherent beams, we need to combine the ``a`` and ``b`` coefficients
before doing the final force calculation. We can either: \* translate
both beams from the coordinate origin to the displaced locations,
combine the beams and then translate the combined beam to the particle
location \* translate each beam to the particle location and combine the
beams before calculating the force

Depending on the size of the beams, the separation and the size of the
particle, these methods will take different amounts of time. If the
particle is significantly smaller than the beam Nmax or the separation
between the two beams, it is most likely faster to combine the beams
after translating the individual beams. If many translations/rotations
are needed, it is most likely faster to create a single beam and apply
the translations to that beam.

Creating a single beam
----------------------

In order to create a single beam, we first need to translate the two
beams from the origin to their displaced locations. In order to be able
to translate a beam multiple times, we need to keep the higher order
multipole terms after the first translation. To calculate the Nmax we
need for this translation, we convert the old Nmax to a radius and add
the displacement before converting back to a Nmax. We then request that
``translateXyz`` produce a beam with the new Nmax.

.. code:: matlab

    % Calculate new Nmax
    Nmax = ott.utils.ka2nmax(ott.utils.nmax2ka(beam.Nmax) ...
        + displacement*T.k_medium);

    % Change the Nmax and create the two beams
    beam1 = beam.translateXyz([-displacement; 0; 0], 'Nmax', Nmax);
    beam2 = beam.translateXyz([displacement; 0; 0], 'Nmax', Nmax);

To combine the beams, we simply add them (which adds the ``a`` and ``b``
coefficients of each beam). We can change the relative phase between the
two beams by simply multiplying a complex phase term by one of the
beams.

.. code:: matlab

    % Add the beams
    nbeam = beam1 + beam2 * phase;

The force can then be calculated from this combined beam:

.. code:: matlab

    % Calculate the force along the x-axis
    fx1 = ott.forcetorque(nbeam, T, 'position', [1;0;0] * x);

Combining after translations/rotations
--------------------------------------

Instead of applying multiple translations to each beam, it is also
possible to apply only a single translation to each beam. There is not
currently any automated method for doing this in the toolbox, the
easiest way is to add the translations and force calculation to a for
loop:

.. code:: matlab

    for ii = 1:length(x)

      % Translate and add the beams
      beam1 = beam.translateXyz([x(ii)+displacement; 0; 0]);
      beam2 = beam.translateXyz([x(ii)-displacement; 0; 0]);
      tbeam = beam1 + beam2 * phase;

      % Scatter the beam and calculate the force
      sbeam = T * tbeam;
      fx2(:, ii) = ott.forcetorque(tbeam, sbeam);
    end

Combining incoherent beams
==========================

For incoherent beams we just need to sum the force from each individual
beam. Similarly to coherent beams, we can apply the translation to the
individual beams or to a combined incoherent beam object.

Operations on individual beams
------------------------------

For incoherent beams, we can use the ``ott.forcetorque`` method to do
the translations and force calculation. Unlike coherent beams, we don't
need to combine the beams after translating the beam.

.. code:: matlab

    fx3 = ott.forcetorque(beam, T, ...
        'position', [1;0;0] * x + [displacement; 0; 0]);
    fx3 = fx3 + ott.forcetorque(beam * phase, T, ...
        'position', [1;0;0] * x - [displacement; 0; 0]);

Applying the same operations on both beams
------------------------------------------

As with coherent beams, combining both beams requires translating the
beam and keeping the higher Nmax terms. We can then combine the beams
into a single ``ott.Bsc`` object and use the ``ott.forcetorque`` method
to apply the translations and calculate the forces. When
``ott.forcetorque`` is called with a ``ott.Bsc`` object containing
multiple beams, it produces a 3-Dimensional matrix with the third
dimension corresponding to the force from each beam in ``ott.Bsc``. To
calculate the incoherent force, we simply need to sum over the third
dimension of this matrix.

.. code:: matlab

    % Calculate new Nmax
    Nmax = ott.utils.ka2nmax(ott.utils.nmax2ka(beam.Nmax) ...
        + displacement*T.k_medium);
      
    % Change the Nmax and create the two beams
    beam1 = beam.translateXyz([-displacement; 0; 0], 'Nmax', Nmax);
    beam2 = beam.translateXyz([displacement; 0; 0], 'Nmax', Nmax);

    beamc = beam1.append(beam2);

    fx4 = ott.forcetorque(beamc, T, 'position', [1;0;0] * x);
    fx4 = sum(fx4, 3);
