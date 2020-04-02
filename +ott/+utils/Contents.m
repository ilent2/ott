% ott.utils utilitiy functions
%
% Files
%   angulargrid                    - makes a angular grid of points over a sphere
%   combined_index                 - translates between (n,m) and combined index
%   hermite                        - calculate hermite polynomial.
%   hgmode                         - calcualtes HG mode amplitude at z = 0
%   incecoefficients               - calculates fourier coefficients for Ince polynomials
%   iseven                         - determines if an integer is even
%   isodd                          - determines if an integer is odd
%   ka2nmax                        - finds a reasonable Nmax to truncate at for given size parameter
%   laguerre                       - associated Laguerre function
%   legendrerow                    - gives the spherical coordinate recursion in m
%   lgmode                         - calculates LG mode amplitude at z = 0
%   matchsize                      - checks that all vector inputs have the same number of rows
%   nmax2ka                        - finds size parameter ka corresponding to Nmax
%   paraxial_beam_waist            - computes the re-normlaised beam waist for high-order
%   paraxial_transformation_matrix - produces paraxial beam mode conversion
%   perpcomponent                  - finds perpendicular (and optionally) parallel
%   rotation_matrix                - calculates rotation matrix using Euler-Rodrigues formula.
%   rtp2xyz                        - coordinate transformation from spherical to cartesian
%   rtpv2xyzv                      - spherical to cartiesn vector field conversion
%   sbesselh                       - spherical hankel function hn(kr) of the first or 
%   sbesselh1                      - spherical hankel function hn(kr) of the first kind,
%   sbesselh2                      - spherical hankel function hn(kr) of the second kind
%   sbesselj                       - spherical bessel function jn(kr)
%   spharm                         - scalar spherical harmonics and angular partial derivatives.
%   threewide                      - creates colum vector with input repeated in 3 columns
%   translate_z                    - calculates translation matricies for translation of
%   vsh                            - calculate vector spherical harmonics
%   vswf                           - vector spherical wavefunctions: M_k, N_k.
%   vswfcart                       - vector spherical harmonics spherical coordinate input,
%   wigner_rotation_matrix         - rotation matrix for rotation of spherical
%   xyz2rtp                        - coordinate transformation from cartesian to spherical
%   xyzv2rtpv                      - cartiesian to spherical vector field conversion
%   rotate_3x3tensor               - rotate a 3x3 tensor or create a 3x3 tensor
%
% Sub-packages
%   +polarizability                - Polarizability calculation methods for DDA
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using this file.
