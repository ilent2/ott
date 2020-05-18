classdef InceGaussian < ott.beam.abstract.Gaussian
% Abstract representation of a Ince-Gaussian beam
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist       -- Beam waist at focus
%   - lmode       -- Azimuthal mode number
%   - porder      -- Paraxial mode order.
%   - ellipticity -- Ellipticity of coordinates
%   - parity      -- Parity of beam ('even' or 'odd')
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
    lmode        % Azimuthal mode number
    porder       % Paraxial mode order.
    ellipticity  % Ellipticity of coordinates
    parity       % Parity of beam ('even' or 'odd')
  end

  methods
    function beam = InceGaussian(varargin)
      % Construct a new Abstract Ince-Gaussian beam
      %
      % Usage
      %   beam = InceGaussian(waist, lmode, porder, parity, ellipticity, ...)
      %
      %   beam = InceGaussian(...)
      %   Construct a beam with `waist = 1`, `lmode = 0`, `porder = 0`,
      %   even parity and `ellipticity = 1`.
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist
      %   - lmode (numeric) -- Azimuthal mode number.
      %   - porder (numeric) -- Paraxial mode number.
      %   - parity (enum) -- Either 'even' or 'odd'.
      %   - ellipticity (numeric) -- Ellipticity of coordinates.
      %
      % For optional parameters, see :class:`Properties`.

      p = inputParser;
      p.addOptional('waist', 1.0, @isnumeric);
      p.addOptional('lmode', 0, @isnumeric);
      p.addOptional('porder', 0, @isnumeric);
      p.addOptional('parity', 'even', @(x) any(strcmpi(x, {'even', 'odd'})));
      p.addOptional('ellipticity', 1, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.abstract.Gaussian(p.Results.waist, unmatched{:});
      beam.parity = p.Results.parity;
      beam.ellipticity = p.Results.ellipticity;
      beam.lmode = p.Results.lmode;
      beam.porder = p.Results.porder;
    end
  end

  methods % Getters/setters
    % lmode        % Azimuthal mode number
    % pmode        % Paraxial mode order.
    % ellipticity  % Ellipticity of coordinates
    % parity       % Parity of beam ('even' or 'odd')

    function beam = set.lmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'lmode must be numeric scalar');
      assert(round(val) == val, 'lmode must be integer');
      beam.lmode = val;
    end

    function beam = set.porder(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'porder must be numeric scalar');
      assert(round(val) == val, 'pmode must be integer');
      beam.porder = val;
    end

    function beam = set.ellipticity(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'ellipticity must be numeric scalar');
      beam.ellipticity = val;
    end

    function beam = set.parity(beam, val)
      assert(any(strcmpi(val, {'even', 'odd'})), ...
          'parity must be ''even'' or ''odd''');
      beam.parity = val;
    end
  end
end

