classdef Ray < ott.beam.properties.Ray ...
    & ott.beam.abstract.Beam
% Specialisation describing abstract geometric ray beams.
% Inherits from :class:`ott.beam.properties.Ray` and :class:`Beam`.
%
% This class is only provided for consistency.  To create a Ray beam,
% it is recommended to create either a :class:`ott.beam.Ray` directly
% or create another beam and cast to a `ott.beam.Ray`, for example
%
% .. code:: matlab
%   beam = ott.beam.abstract.Gaussian(1.0);
%   rays = ott.beam.Ray(beam);
%
% Supported casts
%   - Beam
%   - Ray
%   - PlaneWave
%   - vswf.Bsc
%   - vswf.PlaneWave

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

end
