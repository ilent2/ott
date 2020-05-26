% ott.beam Package containing descriptions of beams.
%
% This package contains descriptions of optical beams.
%
% Beams can be formed into arrays using any of the :class:`Array`,
% :class:`Coherent` or :class:`Incoherent` classes.
% Some classes implicitly support arrays of beams by inheriting from
% :class:`ott.utils.ArrayType`.
%
% Generic beams
%   Beam            -- Base class for all beam implementations
%   Scattered       -- Represents generic scattered beams
%   Array           -- Represents arrays of beams
%   Coherent        -- Pseudonym for coherent arrays
%   Incoherent      -- Pseudonym for incoherent arrays
%   Empty           -- Empty beam array element
%
% Beam specialisations
%   +paraxial       -- Paraxial specialisations sub-package
%   +vswf           -- Vector spherical wave function sub-package
%   PlaneWave       -- Array of plane wave beams
%   Ray             -- Array of geometric rays
%   ScatteredRay    -- Scattered specialisation for Rays
%   GaussianDavis5  -- 5th order approximation for a Gaussian
%   ZeroScattered   -- Scattered beam produced by zero-scattering methods
%   Dipole          -- Radiation field produced by a dipole
%
% Sub-packages
%   +properties     -- Descriptions of beam properties
%   +abstract       -- Abstract representations of beams
%   +medium         -- Descriptions of optical mediums
%   +utils          -- Utilities for describing beams
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

