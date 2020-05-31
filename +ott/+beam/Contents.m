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
%   Array           -- Represents arrays of beams
%
% Beam specialisations
%   +paraxial       -- Paraxial specialisations sub-package
%   +vswf           -- Vector spherical wave function sub-package
%   PlaneWave       -- Array of plane wave beams
%   Ray             -- Array of geometric rays
%   GaussianDavis5  -- 5th order approximation for a Gaussian
%   Dipole          -- Radiation field produced by a dipole
%
% Sub-packages
%   +properties     -- Descriptions of beam properties
%   +abstract       -- Abstract representations of beams
%   +medium         -- Descriptions of optical mediums
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

