classdef CastFarfield < ott.beam.properties.CastFarfield ...
    & ott.beam.abstract.Beam
% Base class with casts for far-field.
% Inherits from :class:`Beam` and :class:`ott.beam.properties.CastFarfield`.
%
% Supported casts
%   - Beam             -- (Inherited) Uses vswf.Bsc
%   - vswf.Bsc         -- (Inherited) Uses vswf.Pointmatch
%   - vswf.Pointmatch  -- Overloaded, Uses vswf.FarfieldPm
%   - vswf.NearfieldPm -- (Inherited) Uses Beam via :class:`CastFarfield`
%   - vswf.FarfieldPm  -- (Inherited)
%   - abstract.InterpFarfield   -- (Inherited)
%   - abstract.InterpNearfield  -- (Inherited) Uses Beam via CastFarfield

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = ott.beam.vswf.Pointmatch(varargin)
      % Cast to vswf.FarfieldPm
      beam = ott.beam.vswf.NearfieldPm(varargin{:});
    end
  end
end
