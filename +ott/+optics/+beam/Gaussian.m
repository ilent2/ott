classdef Gaussian < ott.optics.beam.AbstractBeam ...
    & ott.optics.beam.utils.VariablePower
% Abstract description of a Gaussian beam.
% Inherits from :class:`AbstractBeam` and :class:`utils.VariablePower`.
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
      % For optional parameters, see :class:`BeamProperties`.

      % Filter output power for Variable power constructor
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('power', 1.0);
      p.parse(varargin{:});
      unmatched = [fieldnames(p.Unmatched).', struct2cell(p.Unmatched).'];

      beam = beam@ott.optics.beam.utils.VariablePower(p.Results.power);
      beam = beam@ott.optics.beam.AbstractBeam(unmatched{:});
      beam.waist = waist;
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

