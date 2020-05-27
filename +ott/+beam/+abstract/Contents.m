% ott.beam.abstract Abstract representations of beams
%
% This sub-package contains abstract representations of beams.
% These beams have all the properties of the beam they represent
% but don't implement methods to calculate the fields.  They define
% methods for creating beams of other types (casts).
%
% Abstract beams inherit from ott.beam.properties classes with the
% same names but not :class:`ott.beam.Beam`.  They define casts and
% should not be used as base classes for other beams.
%
% Abstract beams are hetrogeneous arrays and do not inherit from
% :class:`ott.beam.utils.ArrayType`.  Arrays of beams are assumed
% to be coherent unless they are :class:`Array` or :class:`Incoherent`.
% The :class:`Coherent` class is semi-redundant, only useful for
% arrays of sub-arrays of coherent beams.
%
% Generic beams
%   Beam              -- Base class for abstract beam representations.
%   Scattered         -- Represents generic scattered beams
%   Array             -- Represents generic arrays of beams
%   Coherent          -- Represents arrays of coherent beams
%   Incoherent        -- Represents arrays of incoherent beams
%   Empty             -- Empty beam array element
%
% Paraxial beams
%   Gaussian          -- Gaussian beam.
%   HermiteGaussian   -- Hermite-Gaussian beam.
%   LaguerreGaussian  -- Laguerre-Gaussian beam.
%   InceGaussian      -- Ince-Gaussian beam.
%
% PlaneWave-like beams
%   PlaneWave         -- Plane wave beam.
%   Ray               -- Geometric optics ray beam.
%
% Bessel-like beams
%   Bessel            -- Bessel beam
%   Annular           -- Annular beam
%
% Masked beams
%   MaskedFarfield    -- Combination beam with far-field mask
%   MaskedNearfield2d -- Combination beam with 2-D near-field mask
%   MaskedNearfield3d -- Combination beam with 3-D near-field mask
%   MaskedParaxial    -- Combination beam with paraxial mask
%
% Interpolated beams
%   InterpParaxial    -- Interpolated beam for paraxial-field
%   InterpNearfield2d -- Interpolated beam for near-field
%   InterpNearfield3d -- Interpolated beam for near-field
%   InterpFarfield    -- Interpolated beam for far-field
%
% Other beams
%   UniformFarfield   -- Beam with a uniform far-field
%   TopHat            -- Collimated top-hat beam
%   FocussedTopHat    -- Focussed paraxial top-hat beam
%   Dipole            -- Radiation field produced by a dipole
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

