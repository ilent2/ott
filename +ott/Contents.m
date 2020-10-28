% The following was originally notes, but parts should probably become
% the Contents.m file and be added to the documentation.

bsc and tmatrix directories describe the operations themselves

beam and particle directories are the "simple" interface.
Each beam and particle has a position, rotation.  Particles also have
shapes and drag,

We should use the new VswfData class in the beam creation functions, for
example, when we do point matching for Gaussian or BscPmParaxial.
This might significantly speed up beam creation.

We should implement specialisation of Bsc for PlaneWave and Annular.
These beams have more specific translation and rotation implementations.
Should we also have a ott.bsc.scattered specialisation or would this
be better suited as a ott.beam.Scattered class?  Probably ott.beam?

In the previous draft, we also had reduced basis specialisations.
These might be useful too?

We also probably want to implement a GrowOnUse helper or something similar.

T-matrices have no geometry, but particles can (optionally).

Where are we going to put DDA?

