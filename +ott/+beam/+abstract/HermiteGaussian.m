classdef HermiteGaussian < ott.beam.abstract.Gaussian
% Abstract representation of a Hermite-Gaussian beam
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist      -- Beam waist at focus
%   - mmode      -- Hermite mode order
%   - nmode      -- Hermite mode order
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
    mmode       % Hermite mode order
    nmode       % Hermite mode order
  end

  methods
    function beam = HermiteGaussian(waist, mmode, nmode, varargin)
      % Construct a new Abstract Hermite-Gaussian beam
      %
      % Usage
      %   beam = HermiteGaussian(waist, mmode, nmode, ...)
      %
      % Parameters
      %   - waist (numeric) -- Beam waist
      %   - mmode (integer) -- Mode number
      %   - nmode (integer) -- Mode number
      %
      % For optional parameters, see :class:`Properties`.

      beam = beam@ott.beam.abstract.Gaussian(waist, varargin{:});
      beam.mmode = mmode;
      beam.nmode = nmode;
    end
  end

  methods % Getters/setters
    % mmode       % Hermite mode order
    % nmode       % Hermite mode order

    function beam = set.nmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'nmode must be numeric scalar');
      assert(round(val) == val, 'nmode must be integer');
      beam.nmode = val;
    end

    function beam = set.mmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'mmode must be numeric scalar');
      assert(round(val) == val, 'mmode must be integer');
      beam.mmode = val;
    end
  end
end

