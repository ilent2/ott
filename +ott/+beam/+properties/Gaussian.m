classdef Gaussian < ott.beam.properties.Polarisation ...
    & ott.beam.properties.Mapping
% Properties of Gaussian beams
%
% Properties
%   - waist         -- Beam waist radius [m]
%   - mapping       -- Paraxial to far-field beam mapping
%   - polbasis      -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- (2 numeric) Field in theta/phi or x/y directions
%   - truncation_angle -- Truncation angle (hard edge for aperture)
%
% Methods
%   - isGaussian    -- Returns true

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    waist          % Beam waist radius [m]
    truncation_angle % Truncation angle [radians]
  end

  properties (Hidden, SetAccess=protected)
    waistInternal
    truncation_angleInternal
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
      assert(isnumeric(val) && isscalar(val) && val >= 0, ...
        'waist must be positive numeric scalar');
      
      % Warning user about large beam waists
      if val > 1e-3
        warning('ott:beam:properties:Gaussian:large_waist', ...
          'beam waist larger than 1e-3 meters');
      end
      
      beam.waistInternal = val;
      beam.data = [];
    end
    function val = get.waist(beam)
      val = beam.waistInternal;
    end

    function beam = set.truncation_angle(beam, val)
      assert(isnumeric(val) && isscalar(val) && val > 0 && val <= pi, ...
        'truncation_angle must be numeric scalar between [0, pi]');
      beam.truncation_angleInternal = val;
      beam.data = [];
    end
    function val = get.truncation_angle(beam)
      val = beam.truncation_angleInternal;
    end
  end
end

