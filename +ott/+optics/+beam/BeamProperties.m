classdef (Abstract) BeamProperties
% A base class for Beam and AbstractBeam representations
%
% This class defines the common properties and methods to these
% two classes.
%
% Properties
%   - power       -- The power of the beam (may be infinite)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%
% Dependent properties
%   - wavenumber    -- Wave-number of beam in medium
%
% Abstract methods
%   - getBeamPower      -- get method called by dependent property power

  properties
    wavelength     % Wavelength of beam in medium (default: 1.0)
  end

  properties (Dependent)
    power           % The power of the beam (may be infinite)
    wavenumber      % Wave-number of beam in medium
  end

  methods (Hidden)
    function beam = setBeamPower(beam, val)
      % Function to set the beam power (if supported)
      % Override this function if your beam supports this feature
      error('Setting beam power not supported');
    end
  end

  methods
    function val = get.power(beam)
      val = beam.getBeamPower();
    end
    function beam = set.power(beam, val)
      beam = beam.setBeamPower(val);
    end

    function beam = set.wavelength(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'wavelength must be numeric scalar');
      beam.wavelength = val;
    end

    function beam = set.wavenumber(beam, val)
      % Set the wavelength
      assert(isnumeric(val), 'wavenumber must be numeric');
      beam.wavelength = 2*pi./val;
    end
    function val = get.wavenumber(beam)
      val = 2*pi/beam.wavelength;
    end
  end
end
