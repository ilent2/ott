classdef LaguerreGaussian < ott.beam.properties.Gaussian ...
    & ott.beam.properties.Lmode
% Properties of Laguerre-Gaussian beams.
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist         -- Beam waist at focus
%   - pmode         -- Radial Laguerre mode order
%   - lmode         -- (numeric) Orbital angular momentum of the beam
%   - power         -- Beam power
%   - polbasis      -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- (2 numeric) Field in theta/phi or x/y directions
%   - mapping       -- Paraxial to far-field beam mapping
%
% Methods
%   - isGaussian    -- Returns true if lmode == pmode == 0

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    pmode          % Radial Laguerre mode order
  end

  properties (Hidden, SetAccess=protected)
    pmodeInternal
  end

  properties (Abstract)
    data
  end

  methods
    function b = isGaussian(beam)
      b = beam.lmode == 0 && beam.pmode == 0;
    end
  end

  methods % Getters/setters
    function beam = set.pmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'pmode must be numeric scalar');
      assert(round(val) == val && val >= 0, ...
          'pmode must be positive integer');
      beam.pmodeInternal = val;
      beam.data = [];
    end
    function val = get.pmode(beam)
      val = beam.pmodeInternal;
    end
  end
end

