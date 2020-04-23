classdef Gaussian < ott.optics.beam.AbstractBeam
% Description of a Gaussian beam
%
% Units of the fields depend on units used for the properties.
%
% Properties
%   - waist         -- Beam waist radius
%   - permittivity  -- Relative permittivity of medium (default: 1.0)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%   - speed0        -- Speed of light in vacuum (default: 1.0)
%
% Properties (Dependent)
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
    permittivity   % Relative permittivity of medium (default: 1.0)
    wavelength     % Wavelength of beam in medium (default: 1.0)
    speed0         % Speed of light in vacuum (default: 1.0)
  end

  properties (Dependent)
    omega          % Optical frequency of light
    index_medium   % Refractive index in medium
    wavenumber     % Wave-number of beam in medium
    speed          % Speed of light in medium
    wavelength0    % Vacuum wavelength of beam
  end

  methods (Static)

    function [permittivity, wavelength, speed0] = parseInputs(varargin)
      % Helper to parse the inputs to the constructor

      % All values should be numeric scalars
      type_check = @(x) isnumeric(x) & isscalar(x);

      % Setup input arguments
      p = ott.utils.RelatedArgumentParser;
      p.addRequired('permittivity', 1.0, type_check);
      p.addRequired('wavelength', 1.0, type_check);
      p.addRequired('speed0', 1.0, type_check);
      p.addOptional('omega', [], type_check);
      p.addOptional('index_medium', [], type_check);
      p.addOptional('wavenumber', [], type_check);
      p.addOptional('speed', [], type_check);
      p.addOptional('wavelength0', [], type_check);

      % Add rules relating required to optional

      p.addRule('speed0', @(x, y) x .* y, 'wavelength0', 'omega');
      p.addRule('speed0', @(x, y) sqrt(x) .* y, 'permittivity', 'speed');
      p.addRule('speed0', @(x, y, z) sqrt(x) .* y .* z ./ (2*pi), ...
          'permittivity', 'wavelength', 'omega');

      p.addRule('permittivity', @(x) x.^2, 'index_medium');
      p.addRule('permittivity', @(x, y) (x ./ y).^2, 'speed0', 'speed');
      p.addRule('permittivity', @(x, y) (x ./ y).^2, ...
          'wavelength0', 'wavelength');
      p.addRule('permittivity', @(x, y, z) (x ./ y ./ z .* (2*pi)).^2, ...
          'speed0', 'wavelength', 'omega');

      p.addRule('wavelength', @(x) 2.*pi./x, 'wavenumber');
      p.addRule('wavelength', @(x, y) x ./ sqrt(y), ...
          'wavelength0', 'permittivity');
      p.addRule('wavelength', @(x, y) x ./ y .* (2*pi), 'speed', 'omega');
      p.addRule('wavelength', @(x, y, z) x ./ sqrt(y) ./ z .* (2*pi), ...
          'speed0', 'permittivity', 'omega');

      % Parse inputs
      p.parse(varargin{:});

      % Assign variables to output
      speed0 = p.RequiredResults.speed0;
      permittivity = p.RequiredResults.permittivity;
      wavelength = p.RequiredResults.wavelength;
    end
  end

  methods
    function beam = Gaussian(waist, varargin)
      % Construct a new Gaussian beam description
      %
      % Usage
      %   beam = Gaussian(waist, ...)
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [L]
      %
      % Optional named arguments
      %   - permittivity (numeric) -- Relative permittivity of medium
      %   - wavelength (numeric) -- Wavelength in medium [L]
      %   - speed0 (numeric) -- Speed of light in vacuum [L/T]
      %
      %   - omega (numeric) -- Optical frequency [2*pi/T]
      %   - index_medium (numeric) -- Refractive index in medium
      %   - wavenumber (numeric) -- Wave-number in medium [2*pi/L]
      %   - speed (numeric) -- Speed of light in medium [L/T]
      %   - wavelength0 (numeric) -- Wavelength in medium [L]

      % Store waist
      beam.waist = waist;

      % Parse optional inputs
      [beam.permittivity, beam.wavelength, beam.speed0] = ...
          beam.parseInputs(varargin{:});
    end
  end

  methods
    function beam = set.waist(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'waist must be numeric scalar');
      beam.waist = val;
    end
    function beam = set.wavelength(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'wavelength must be numeric scalar');
      beam.wavelength = val;
    end
    function beam = set.speed0(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'speed0 must be numeric scalar');
      beam.speed0 = val;
    end
    function beam = set.permittivity(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'permittivity must be numeric scalar');
      beam.permittivity = val;
    end

    function val = get.omega(beam)
      wavelength0 = beam.wavelength .* beam.index_medium;
      val = 2*pi*beam.speed0./wavelength0;
    end
    function val = get.index_medium(beam)
      val = sqrt(beam.permittivity);
    end
    function val = get.wavenumber(beam)
      val = 2*pi/beam.wavelength;
    end
    function val = get.speed(beam)
      val = beam.speed0 ./ beam.index_medium;
    end
    function val = get.wavelength0(beam)
      val = beam.wavelength ./ beam.index_medium;
    end
  end
end

