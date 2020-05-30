classdef HomogeneousRelative < ott.scat.Particle
% Declares a index_relative property for homogeneous particles.
%
% Properties
%   - index_relative      -- Relative refractive index

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    index_relative      % Relative refractive index of particle
  end

  methods
    function particle = HomogeneousRelative(index_relative)
      % Construct and assign the relative refractive index
      %
      % Usage
      %   particle = HomogeneousRelative(index_relative)

      particle.index_relative = index_relative;
    end
  end

  methods % Getters/setters
    function particle = set.index_relative(particle, val)
      assert(isnumeric(val) && isscalar(val), ...
          'index_relative must be numeric scalar');
      particle.index_relative = val;
    end
  end
end
