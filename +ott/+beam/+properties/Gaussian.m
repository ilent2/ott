classdef Gaussian < ott.beam.properties.Polarisation ...
    & ott.beam.properties.Mapping
% Properties of Gaussian beams
%
% Properties
%   - waist         -- Beam waist radius [m]
%   - power         -- Beam power [W]
%   - mapping       -- Paraxial to far-field beam mapping
%   - polbasis      -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- (2 numeric) Field in theta/phi or x/y directions
%
% Methods
%   - isGaussian    -- Returns true

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    waist          % Beam waist radius [m]
    power          % Beam power [W]
  end

  properties (Hidden, SetAccess=protected)
    waistInternal
    powerInternal
  end

  properties (Abstract)
    data
  end

  methods
    function b = isGaussian(~)
      b = true;
    end
  end

  methods % Getters/setters
    function beam = set.waist(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'waist must be numeric scalar');
      beam.waistInternal = val;
      beam.data = [];
    end
    function val = get.waist(beam)
      val = beam.waistInternal;
    end

    function beam = set.power(beam, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'power must be positive numeric scalar');
      beam.powerInternal = val;
      beam.data = [];
    end
    function val = get.power(beam)
      val = beam.powerInternal;
    end
  end
end

