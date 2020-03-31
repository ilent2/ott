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
%
% Files for DDA
%   E_inc_vswf            - Calculate the incident field for a single vswf mode
%   calc_Aj               - calculates a 3 X 3N block comprising N number of 3 X 3 Green's tensors  
%   interaction_A         - Calculate the interaction matrix
%   polarizability_CM     - Clausius-Mossoti Polarizability 
%   polarizability_FCD    - Filtered coupled dipole polarizability
%   polarizability_LDR    - Lattice dispersion relation polarizablity
%   rotate_polarizability - Apply a set of rotations to the unit-polarisability
%   col3to1               - Reshape 3 column vector to 1 column vector format
%   col1to3               - Reshape 1 column vector to 3 column vector format

