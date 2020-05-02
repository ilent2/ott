classdef Dielectric < ott.optics.beam.dipole.Dipole
% A dielectric dipole instance
%
% The dipole polarizability is calculated from the particle relative
% refractive index.
%
% Properties
%   - index_relative      -- Relative refractive index
%   - index_medium        -- Refractive index of medium (default: 1.0)
%
% Dependent properties
%   - index_particle      -- Refractive index of particle

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    index_relative       % Relative refractive index
    index_medium         % Refractive index of medium (default: 1.0)
  end

  properties (Dependent)
    index_particle       % Refractive index of particle
  end

  methods
    function dipole = Dielectric()
      % TODO
      error('not implemented');
    end
  end

  methods
    function dipole = set.index_relative(dipole, val)
      assert(isnumeric(val) && isscalar(val), ...
        'index_relative must be numeric scalar');
      dipole.index_relative = val;
    end
    function dipole = set.index_medium(dipole, val)
      assert(isnumeric(val) && isscalar(val), ...
        'index_medium must be numeric scalar');
      dipole.index_medium = val;
    end

    function val = get.index_particle(dipole)
      % Calculate index medium from relative and particle index
      % index_relative = index_particle / index_medium
      val = dipole.index_relative .* dipole.index_medium;
    end
  end
end
