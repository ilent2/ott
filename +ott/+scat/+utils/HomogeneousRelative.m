classdef HomogeneousRelative
% Declares a index_relative property for homogeneous particles.
%
% Properties
%   - index_relative      -- Relative refractive index
%
% No class constructor is implemented.  Property is not initialised.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    index_relative      % Relative refractive index of particle
  end

  methods % Getters/setters
    function particle = set.index_relative(particle, val)
      assert(isnumeric(val) && isscalar(val), ...
          'index_relative must be numeric scalar');
      particle.index_relative = val;
    end
  end
end
