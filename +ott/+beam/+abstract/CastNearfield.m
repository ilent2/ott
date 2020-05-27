classdef CastNearfield < ott.beam.properties.CastNearfield ...
    & ott.beam.abstract.Beam
% Base class with casts for near-field.
% Inherits from :class:`Beam` and :class:`ott.beam.properties.CastNearfield`.
%
% Supported casts
%   - Beam             -- (Inherited) Uses vswf.Bsc
%   - vswf.Bsc         -- (Inherited) Uses vswf.Pointmatch
%   - vswf.Pointmatch  -- Overloaded, Uses vswf.NearfieldPm
%   - vswf.NearfieldPm -- (Inherited)
%   - vswf.FarfieldPm  -- (Inherited) Uses Beam via :class:`CastNearfield`
%   - abstract.InterpNearfield  -- (Inherited)
%   - abstract.InterpFarfield   -- (Inherited) Uses Beam via CastNearfield

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = ott.beam.vswf.Pointmatch(varargin)
      % Cast to vswf.NearfieldPm
      beam = ott.beam.vswf.NearfieldPm(varargin{:});
    end
  end
end
