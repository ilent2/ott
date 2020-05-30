% ott.beam.properties Descriptions of beam properties
%
% This sub-package contains classes describing properties of beams
% which can be re-used by various beam implementations.  These classes
% don't inherit from any special classes such as hetrogeneous array,
% making them compatible with both beam and abstract beam representations.
%
% Generic classes
%   Beam              -- Base class for beam properties.
%   Material          -- Base class for beams with materials.
%   FarfieldMapping   -- Adds properties and methods for far-field mapping
%   Empty             -- Empty beam array element
%
% Power property
%   InfinitePower     -- Infinite power property
%   VariablePower     -- Variable power property
%   ZeroPower         -- Zero power property
%
% Cast methods
%   CastBoth          -- Cast beam for field calculation
%   CastFarfield      -- Cast beam for far-field calculation
%   CastNearfield     -- Cast beam for near-field calculation
%
% Scattered classes
%   ScatteredInterface  -- Interface for scattered properties
%   ZeroScattered       -- Zero-scattered properties
%   Scattered           -- Generic scattered beam properties
%
% Dipole classes
%   Dipole
%   DipoleArray         -- Specialisation inheriting from ArrayType
%
% Plane wave classes
%   PlaneWave           -- Plane wave properties
%   PlaneWaveArray      -- Specialisation inheriting from ArrayType
%
% Array types
%   ArrayType           -- Base class for arrays of beams (non-abstract)
%   AbstractArray       -- Base class for abstract arrays of beams
%   CoherentArrayType   -- Specialisation for coherent array types
%   IncoherentArrayType -- Specialisation for incoherent array types
%   ArrayArrayType      -- Specialisation for generic array types
%   AnyArrayType        -- Specialisation for variable array types
%
% Paraxial beams
%   Gaussian          -- Gaussian beam.
%   HermiteGaussian   -- Hermite-Gaussian beam.
%   LaguerreGaussian  -- Laguerre-Gaussian beam.
%   InceGaussian      -- Ince-Gaussian beam.
%
% Other beams
%   Bessel            -- Bessel beam
%   Annular           -- Annular beam
%   TopHat            -- Collimated top-hat beam
%   FocussedTopHat    -- Focussed paraxial top-hat beam
%
% Masked beams
%   MaskedBeam        -- Base class for masked beams
%   MaskedParaxial    -- Properties for masked paraxial beams
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

