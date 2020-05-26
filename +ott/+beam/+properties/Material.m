classdef Material
% Defines material properties for beams.
%
% Properties
%   - omega       -- Optical frequency of beam.
%   - medium      -- Properties of optical medium.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    omega        % Optical frequency of beam.
    medium       % Properties of optical medium.
  end

  methods (Static)
    function args = likeProperties(other, args)
      if isa(other, 'ott.beam.properties.Material')
        args = ott.utils.addDefaultParameter('omega', other.omega, args);
        args = ott.utils.addDefaultParameter('medium', other.medium, args);
      end
    end
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
