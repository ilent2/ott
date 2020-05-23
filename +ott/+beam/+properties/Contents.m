% ott.beam.properties Descriptions of beam properties
%
% This sub-package contains classes describing properties of beams
% which can be re-used by various beam implementations.  These classes
% don't inherit from any special classes such as hetrogeneous array,
% making them compatible with both beam and abstract beam representations.
%
% Generic classes
%   Beam              -- Base class for beam properties.
%   Scattered         -- Represents generic scattered beams
%   Array             -- Represents arrays of beams
%   Coherent          -- Specialisation for coherent arrays of beams
%   Incoherent        -- Specialisation for incoherent arrays of beams
%   Empty             -- Empty beam array element
%
% Paraxial beams
%   Gaussian          -- Gaussian beam.
%   HermiteGaussian   -- Hermite-Gaussian beam.
%   LaguerreGaussian  -- Laguerre-Gaussian beam.
%   InceGaussian      -- Ince-Gaussian beam.
%
% Plane wave beams
%   PlaneWave         -- Plane wave beam.
%   Ray               -- Geometric optics ray beam.
%
% Other beams
%   Bessel            -- Bessel beam
%   Annular           -- Annular beam
%   TopHat            -- Collimated top-hat beam
%   FocussedTopHat    -- Focussed paraxial top-hat beam
%   Dipole            -- Radiation field produced by a dipole
%   FarfieldMasked    -- Combination beam with far-field mask
%   NearfieldMasked   -- Combination beam with near-field mask
%   ParaxialMasked    -- Combination beam with paraxial mask
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

