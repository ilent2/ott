classdef CastBoth < ott.beam.properties.CastNearfield ...
    & ott.beam.properties.CastFarfield ...
    & ott.beam.abstract.Beam
% Base class for beams with casts for both near/far field calculation
% Inherits from :class:`ott.beam.properties.CastNearfield`
% and :class:`ott.beam.properties.CastFarfield` and :class:`Beam`.
%
% Supported casts
%   - Beam             -- (Inherited)
%   - vswf.Bsc         -- (Inherited)
%   - vswf.Pointmatch  -- (Inherited)
%   - vswf.NearfieldPm -- (Inherited) Uses Beam via :class:`CastNearfield
%   - vswf.FarfieldPm  -- (Inherited) Uses Beam via :class:`CastNearfield

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

end
