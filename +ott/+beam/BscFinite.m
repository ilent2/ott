classdef BscFinite < ott.beam.BscBeam
% Bsc with stored finite representation of the beam.
% Inherits from :class:`BscBeam`.
%
% This class describes beams which can be represented using a finite VSWF
% expansion. The class stores a :class:`Bsc` instance internally.
% Fields at any other location can be found by applying a translation
% matrix to the beam data.  A translated beam will typically not be finite
% unless the VSWF expansion includes the original VSWF region.
%
% Methods
%   - getData       -- Get data for specific Nmax
%   - recalculate   -- Recalculate the beam data

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

end
