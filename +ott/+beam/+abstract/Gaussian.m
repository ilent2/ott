classdef Gaussian < ott.beam.abstract.Beam ...
    & ott.beam.utils.VariablePower
% Abstract description of a Gaussian beam.
% Inherits from :class:`abstract.Beam` and :class:`utils.VariablePower`.
%
% Units of the fields depend on units used for the properties.
%
% Properties
%   - waist         -- Beam waist radius
%
% Inherited properties
%   - permittivity  -- Relative permittivity of medium (default: 1.0)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%   - speed0        -- Speed of light in vacuum (default: 1.0)
%   - omega         -- Optical frequency of light
%   - index_medium  -- Refractive index in medium
%   - wavenumber    -- Wave-number of beam in medium
%   - speed         -- Speed of light in medium
%   - wavelength0   -- Vacuum wavelength of beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    waist          % Beam waist radius
  end

  methods
    function beam = Gaussian(waist, varargin)
      % Construct a new Abstract Gaussian beam
      %
      % Usage
      %   beam = Gaussian(waist, ...)
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [L]
      %
      % Optional named arguments
      %   - power (numeric) -- beam power.  Default: ``1.0``.
      %
      % For optional parameters, see :class:`Properties`.

      % Filter output power for Variable power constructor
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('power', 1.0);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.utils.VariablePower(p.Results.power);
      beam = beam@ott.beam.abstract.Beam(unmatched{:});
      beam.waist = waist;
    end
    
    function beam = ott.beam.Beam(oldbeam, varargin)
      % Cast the beam to a ott.beam.Beam object
      %
      % The default beam is a ott.beam.GaussianDavis5, since it
      % is a good compromise between speed and accuracy.
      
      beam = ott.beam.GaussianDavis5(oldbeam, varargin{:});
    end
    
    function beam = ott.beam.GaussianDavis5(oldbeam, varargin)
      % Cast beam to a GaussianDavis5
      
      beam = ott.beam.GaussianDavis5(oldbeam.waist, ...
        'omega', oldbeam.omega, 'medium', oldbeam.medium, ...
        'power', oldbeam.power, 'position', oldbeam.position, ...
        'rotation', oldbeam.rotation, varargin{:});
    end
  end

  methods % Getters/setters
    function beam = set.waist(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'waist must be numeric scalar');
      beam.waist = val;
    end
  end
end

