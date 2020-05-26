classdef (Abstract) MaterialBeam < ott.beam.properties.Beam
% Base class for beams defining material properties.
% Inherits from :class:`Beam` and :class:`ott.beam.properties.MaterialBeam`.
%
% Properties
%   - omega       -- Optical frequency of beam.
%   - medium      -- Properties of optical medium.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  Properties
    omega        % Optical frequency of beam.
    medium       % Properties of optical medium.
  end

  methods % Getters/setters
    function beam = set.omega(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'omega must be numeric scalar');
      beam.omega = val;
    end

    function beam = set.medium(beam, val)
      assert(isa(val, 'ott.beam.medium.Medium'), ...
          'medium must be a medium.Medium');
      beam.medium = val;
    end
  end
end
