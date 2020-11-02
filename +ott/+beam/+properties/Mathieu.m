classdef Mathieu < ott.beam.properties.Annular & ott.beam.properties.Parity
% Properties for Mathieu beams
%
% Properties
%   - theta       -- Annular angle [radians]
%   - morder      -- Mathieu beam mode number
%   - ellipticity -- Ellipticity of Mathieu beam
%   - parity      -- Parity of beam ('even' or 'odd')

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    morder       % Mathieu beam mode number
    ellipticity  % Ellipticity of Mathieu beam
  end

  properties (Hidden, SetAccess=protected)
    morderInternal
    ellipticityInternal
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.morder(beam, val)
      assert(isnumeric(val) && isscalar(val) && round(val) == val, ...
          'morder must be numeric integer');
      beam.morderInternal = val;
      beam.data = [];
    end
    function val = get.morder(beam)
      val = beam.morderInternal;
    end

    function beam = set.ellipticity(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'ellipticity must be numeric scalar');
      beam.ellipticityInternal = val;
      beam.data = [];
    end
    function val = get.ellipticity(beam)
      val = beam.ellipticityInternal;
    end
  end
end

