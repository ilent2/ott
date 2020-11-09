% ott.beam Package providing optical tweezers beams.
%
% The beam package provides representations of optical tweezers beams.
% Classes in this package can be used to construct VSWF representations
% of typical beams used in OT.  These classes provide a slightly simplified
% interface to directly using the ott.bsc classes (along with additional
% features for smart Nmax selection).

% Base classes
%   Beam              - Base class for beam representations
%   Empty             - A beam with no fields
%   Array             - Arrays of incoherent/coherent beams
%   Scattered         - Description of beams scattered by particles
%   ArrayType         - Base class for beam array types
%   BscBeam           - Base class for BSC beams
%   BscFinite         - Base class for finite BSC beams
%   BscInfinite       - Base class for infinite Bsc beams
%   BscWBessel        - Base class for Annular/Bessel like beams
%
% Beam types
%   Gaussian          - Gaussian beam
%   HermiteGaussian   - Hermite--Gaussian beam
%   LaguerreGaussian  - Laguerre-Gaussian beam
%   InceGaussian      - Ince-Gaussian beam
%   PlaneWave         - Plane wave
%   Bessel            - Bessel beam
%   Webber            - Webber beam
%   Mathieu           - Mathieu beam
%   Annular           - Annular beam
%   PmParaxial        - Paraxial point matched beam
%
% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

