% ott.beam.medium Describes material properties of beam media
%
% The toolbox uses four different types of classes to declare optical
% media: mediums, vacuums, materials and relative materials.
% Relative materials define the ratio between mediums/materials and
% are mainly used by the scattering methods.
% Materials declare dimensionless material properties as ratios of the
% vacuum medium (i.e. materials are declared with relative properties).
% Vacuums declare the relative units to use for beam properties, units
% must have SI dimensions (but can be scaled to have different units).
% Mediums combine a vacuum, material and an optical frequency and are
% used by the beam objects.
%
% Mediums
%   Medium    -- Defines an optical medium (material + frequency + units)
%
% Units (vacuums)
%   Vacuum    -- Declares the vacuum medium properties (units)
%   Vacuum.Unitary    -- Vacuum with unitary permittivity/permeability
%   Vacuum.BaseSi     -- Vacuum using base SI units
%
% Materials
%   Material    -- Base class for material definitions
%   Dielectric  -- Material with unitary permeability
%   Arbitrary   -- Arbitrary magnetic/electronic material
%   Generic     -- Generic material definitions
%   Relative    -- Describes a ratio between two materials
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

