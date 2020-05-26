classdef Dipole < ott.beam.properties.MaterialBeam
% Properties of dipoles.
% Inherits from :class:`ott.beam.properties.Beam`.
%
% Properties
%   - position      -- Position of the dipole(s)
%   - rotation      -- Orientation of the dipole(s)
%
% Abstract properties
%   - power         -- Power of the dipole radiated field
%   - polarization  -- Polarization of the dipole(s)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract)
    polarization  % Polarization of the dipole(s)
  end
end
