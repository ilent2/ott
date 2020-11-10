classdef InceGaussian < ott.beam.properties.Gaussian ...
    & ott.beam.properties.Parity
% Properties of Ince-Gaussian beams
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist         -- Beam waist at focus
%   - lmode         -- Azimuthal mode number
%   - porder        -- Paraxial mode order.
%   - ellipticity   -- Ellipticity of coordinates
%   - parity        -- Parity of beam ('even' or 'odd')
%   - power         -- Beam power
%   - polbasis      -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- (2 numeric) Field in theta/phi or x/y directions
%   - mapping       -- Paraxial to far-field beam mapping
%
% Methods
%   - isGaussian    -- Returns true if mode is Gaussian

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    lmode          % Azimuthal mode number
    porder         % Paraxial mode order.
    ellipticity    % Ellipticity of coordinates
  end

  properties (Hidden, SetAccess=protected)
    lmodeInternal
    porderInternal
    ellipticityInternal
  end

  properties (Abstract)
    data
  end

  methods
    function b = isGaussian(beam)
      if beam.lmode ~= 0
        b = false;
      else
        % TODO: Implement this
        error('not yet implemented');
      end
    end
  end

  methods % Getters/setters
    function beam = set.lmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'lmode must be numeric scalar');
      assert(round(val) == val, 'lmode must be integer');
      beam.lmodeInternal = val;
      beam.data = [];
    end
    function val = get.lmode(beam, val)
      val = beam.lmodeInternal;
    end

    function beam = set.porder(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'porder must be numeric scalar');
      assert(round(val) == val, 'pmode must be integer');
      beam.porderInternal = val;
      beam.data = [];
    end
    function val = get.porder(beam)
      val = beam.porderInternal;
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

