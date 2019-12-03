
* Consider alternatives to ott.warning

* Clean up documentation in all files, make sure it is all similar

* More thorough unit tests

* Add anisotropic and inhomogeneous DDA functionality
* Add far-field DDA/T-matrix functionality

* Add DDA and SMARTIES as a supported method in Tmatrix.simple

* Merge forcetorque into Bsc and Tmatrix classes
  This will make it more efficient, since we don't need to have
  weird cases for torque/force/spin, we just ask for what we want.

* Check OTTv2 directory and phd2018/phd2019 directories for other
  things we want to include in this version
