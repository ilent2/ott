classdef EccentricSpheresNn < ott.drag.Stokes
% Calculate drag on an eccentric sphere using Gibson's NN approach.
%
% Uses the NN from
%
%   Lachlan J. Gibson, et al. Phys. Rev. E 99, 043304
%   https://doi.org/10.1103/PhysRevE.99.043304
%
% Properties
%   - innerRadius      -- Radius of inner sphere
%   - outerRadius      -- Radius of outer sphere
%   - separation       -- Minimum separation between spheres

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    innerRadius       % Radius of inner sphere
    outerRadius       % Radius of outer sphere
    separation        % Minimum separation between spheres
  end

  methods
  end
end

