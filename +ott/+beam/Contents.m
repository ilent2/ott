% ott.beam Package containing descriptions of beams.
%
% This package contains descriptions of optical beams.
%
% Generic beams
%   Beam            -- Base class for all beam implementations
%   Scattered       -- Represents generic scattered beams
%   Array           -- Represents arrays of beams
%   Coherent        -- Specialisation for coherent arrays of beams
%   Incoherent      -- Specialisation for incoherent arrays of beams
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

