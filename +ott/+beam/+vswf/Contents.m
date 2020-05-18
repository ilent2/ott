% ott.beam.vswf Vector spherical wave function representations of beams.
%
% Descriptions of beams in a vector spherical wave function (VSWF) basis.
% The Bsc class encapsulates VSWF coefficients, other classes provide
% methods for generating Bsc of various beams.
%
% Generic classes
%   Bsc               - Class representing VSWF beam shape coefficients
%   Scattered         - A Bsc instance describing scattered beams
%   Pointmatch        - A Bsc instance with methods for point-matching
%   FarfieldPm        - Far-field point matching base class
%   NearfieldPm       - Near-field point matching base class
%
% Reduced basis specialisations
%   BesselBasis       - Bsc using a Bessel function basis set
%   PlaneBasis        - Bsc using a Plane wave basis set
%   LgParaxialBasis   - Bsc using a Paraxial LG basis set
%
% Beam specialisations
%   Gaussian          - Bsc specialisation for Gaussian beams
%   LaguerreGaussian  - Bsc specialisation for Laguerre-Gaussian beams
%   HermiteGaussian   - Bsc specialisation for Hermite-Gaussian beams
%   InceGaussian      - Bsc specialisation for Ince-Gaussian beams
%   PlaneWave         - Bsc specialisation for plane wave beams
%   Bessel            - Bsc specialisation for Bessel beams
%   Annular           - Bsc specialisation for annular beams
%   PlaneWave         - Bsc specialisation for plane-wave beams
%
% Utility classes
%   GrowOnUse         - Base class for beams whose Nmax expands when used
%   BscScalar         - Helper for non-Bsc-array beams
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

