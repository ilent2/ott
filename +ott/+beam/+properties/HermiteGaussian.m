classdef HermiteGaussian < ott.beam.properties.Gaussian
% Properties of Hermite-Gaussian beams.
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist         -- Beam waist at focus
%   - mmode         -- Hermite mode order
%   - nmode         -- Hermite mode order
%   - power         -- Beam power
%   - polbasis      -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- (2 numeric) Field in theta/phi or x/y directions
%   - mapping       -- Paraxial to far-field beam mapping
%
% Methods
%   - isGaussian    -- Returns true if mmode == nmode == 0

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    mmode          % Hermite mode order
    nmode          % Hermite mode order
  end

  properties (SetAccess=protected, Hidden)
    mmodeInternal
    nmodeInternal
  end

  properties (Abstract)
    data
  end

  methods
    function b = isGaussian(beam)
      b = beam.mmode == 0 && beam.nmode == 0;
    end
  end

  methods % Getters/setters
    function beam = set.nmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'nmode must be numeric scalar');
      assert(round(val) == val, 'nmode must be integer');
      beam.nmodeInternal = val;
      beam.data = [];
    end
    function val = get.nmode(beam)
      val = beam.nmodeInternal;
    end

    function beam = set.mmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'mmode must be numeric scalar');
      assert(round(val) == val, 'mmode must be integer');
      beam.mmodeInternal = val;
      beam.data = [];
    end
    function val = get.mmode(beam)
      val = beam.mmodeInternal;
    end
  end
end
