% ott.utils -- Support functions and classes for OTT
%
% This directory contains various functions and classes used by components
% of the optical tweezers toolbox.  This includes the vector spherical
% wave function methods, functions for working with different coordinate
% systems, and various other helpful functions.
%
% Classes
%   BesselData      - Stores pre-computed Bessel data for field calculations
%   SpharmData      - Stores pre-computed Spharm data for field calculations
%   VswfData        - Stores pre-computed VSWF data for field calculations
%   FieldVector     - Base class for classes encapsulating field vector data
%   FieldVectorCart - Cartesian field vector instance
%   FieldVectorSph  - Spherical field vector instance
%   RotateHelper    - Helper class providing many short-hand rotation methods
%   RotationPositionEntity - Rotation/position interface for OTT entities
%   RotationPositionProp   - Adds rotation and position properties to a class
%   TranslateHelper        - A helper class providing translation methods
%
% Functions
%   angulargrid        - makes a angular grid of points over a sphere
%   arrayMaxSize       - Get maximum array size that fits in physical memory
%   cart2sph_mat       - Cartesian to spherical coordinate conversion matrix
%   combined_index     - translates between (n,m) and combined index
%   crossProductMatrix - Calculate the cross product matrix
%   incecoefficients   - calculates fourier coefficients for Ince polynomials
%   legendrerow        - gives the spherical coordinate recursion in m
%   inpolyhedron       - Tests if points are inside a 3D triangulated surface
%   ka2nmax            - Calculates reasonable VSWF/Nmax for size parameter
%   laguerre           - associated Laguerre function
%   matchsize          - Checks/matches number of rows in all inputs
%   na2angle           - Converts from numerical aperture to half-angle
%   nargoutCheck       - Checks number of output arguments for class type
%   paraxial2rtp       - Convert from paraxial to spherical coordinates.
%   paraxial_transformation_matrix - produces paraxial beam mode conversion
%   rotx               - Build 3x3 rotation matrix for rotation about x axis
%   rotz               - Build 3x3 rotation matrix for rotation about z axis
%   rtp2xyz            - Coordinate transformation from spherical to Cartesian
%   rtpv2xyzv          - Spherical to Cartesian vector field conversion
%   sbesselh1          - spherical hankel function hn(kr) of the first kind,
%   rtpFarfield        - Adds/removes radius in far-field spherical coordinates
%   nmax2ka            - Finds size parameter ka corresponding to VSWF/Nmax
%   rotate_3x3tensor   - Apply a set of rotations to a 3x3 tensor.
%   roty               - Build 3x3 rotation matrix for rotation about y axis
%   sbesselh2          - spherical hankel function hn(kr) of the second kind
%   sbesselj           - spherical bessel function jn(kr)
%   sph2cart_mat       - Spherical to Cartesian coordinate conversion matrix
%   translate_z        - Calculates translation matrices for VSWF along z-axis
%   unmatchedArgs      - Collect inputParser unmatched arguments into cell
%   vsh                - calculate vector spherical harmonics
%   vswf               - vector spherical wavefunctions: M_k, N_k.
%   vswfcart           - vector spherical harmonics spherical coordinate input,
%   spharm             - scalar spherical harmonics and angular derivatives
%   wigner_rotation_matrix - rotation matrix for rotation of spherical
%   xyz2rtp            - Cartesian to spherical coordinate transform
%   xyzv2rtpv          - Cartesian to spherical vector field conversion
%
% Copyright 2018-2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

