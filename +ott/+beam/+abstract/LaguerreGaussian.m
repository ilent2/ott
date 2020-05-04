classdef LaguerreGaussian < ott.optics.beam.Gaussian
% Abstract representation of a Laguerre-Gaussian beam.
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist      -- Beam waist at focus
%   - lmode      -- Azimuthal Laguerre mode order
%   - pmode      -- Radial Laguerre mode order
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
    lmode       % Azimuthal Laguerre mode order
    pmode       % Radial Laguerre mode order
  end

  methods
    function beam = LaguerreGaussian(waist, lmode, pmode, varargin)
      % Construct a new Abstract Laguerre Gaussian beam
      %
      % Usage
      %   beam = LaguerreGaussian(waist, lmode, pmode, ...)
      %
      % Parameters
      %   - waist (numeric)     -- Beam waist
      %   - lmode (integer)     -- Azimuthal LG mode
      %   - pmode (integer > 0) -- Radial LG mode
      %
      % For optional parameters, see :class:`Properties`.

      beam = beam@ott.optics.beam.Gaussian(waist, varargin{:});
      beam.lmode = lmode;
      beam.pmode = pmode;
    end
  end

  methods % Getters/setters
    % waist       % Beam waist at focus
    % lmode       % Azimuthal Laguerre mode order
    % pmode       % Radial Laguerre mode order

    function beam = set.lmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'lmode must be numeric scalar');
      assert(round(val) == val, ...
          'lmode must be integer');
      beam.lmode = val;
    end

    function beam = set.pmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'pmode must be numeric scalar');
      assert(round(val) == val && val >= 0, ...
          'pmode must be positive integer');
      beam.pmode = val;
    end
  end
end

