classdef Medium
% Specification of an optical medium (material + frequency + units)
%
% Properties
%   - material      -- Description of the material
%   - frequency     -- Optical (angular) frequency
%   - vacuum        -- Vacuum (describes units)
%
% Static properties (methods)
%   - DefaultFrequency  -- Default optical frequency
%   - DefaultVacuum     -- Default vacuum
%
% Dependent properties
%   - permittivity  -- Permittivity of medium
%   - permeability  -- Permeability of medium
%   - wavelength    -- Wavelength in medium
%   - wavelength0   -- Wavelength in vacuum
%   - speed         -- Wave speed in medium
%   - index         -- Refractive index of medium
%   - impedance     -- Impedance of medium
%
% Methods
%   - rdivide       -- Create a relative material
%
% See also :class:`ott.beam.medium.Vacuum`, :class:`ott.beam.medium.Material`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    material       % Description of the material
    frequency      % Optical (angular) frequency
    vacuum         % Vacuum (describes units)
  end

  properties (Dependent)
    permittivity
    permeability
    wavelength      % Wavelength in medium
    wavelength0     % Wavelength in vacuum
    speed           % Speed in medium
    index
    impedance
  end

  methods (Static)
    function out = DefaultFrequency(data)
      % Static variable for default optical frequency
      %
      % Usage
      %   frequency = Medium.DefaultFrequency
      %
      %   Medium.DefaultFrequency(value)
      %   -- Change the default optical frequency used in the constructor.
      %
      % Initial value for default frequency is ``2*pi``.

      persistent Frequency;

      if nargin == 0
        if isempty(Frequency)
          Frequency = 2*pi;
        end
        out = Frequency;
      else
        assert(isnumeric(data) && isscalar(data), ...
            'value must be numeric scalar');
        Frequency = data;
      end
    end

    function out = DefaultVacuum(data)
      % Get/set default vacuum static variable
      %
      % Usage
      %   vacuum = Medium.DefaultVacuum
      %
      %   Medium.DefaultVacuum(value)
      %   -- Change the default vacuum used in the constructor.
      %
      % Default initial vacuum is ``ott.beam.medium.Vacuum.Unitary``.

      persistent Vacuum;

      if nargin == 0
        if isempty(Vacuum)
          Vacuum = ott.beam.medium.Vacuum.Unitary;
        end
        out = Vacuum;
      else
        assert(isa(data, 'ott.beam.medium.Vacuum'), ...
            'value must be a ott.beam.medium.Vacuum');
        Vacuum = data;
      end
    end

    function med = FromWavelength(varargin)
      % Construct medium using wavelength instead of frequency
      %
      % Usage
      %   medium = Medium.FromWavelength(material, wavelength, vacuum)
      %
      % Parameters
      %   - material (Material) -- Material properties of medium.
      %
      %   - wavelength (numeric) -- Wavelength in medium.
      %
      %   - vacuum (Vacuum) -- Vacuum (defines units).
      %     Default: ``Medium.DefaultVacuum()``.

      p = inputParser;
      p.addOptional('material', [], @(x) isa(x, 'ott.beam.medium.Material'));
      p.addOptional('wavelength', [], @(x) ~isempty(x) && isnumeric(x));
      p.addOptional('vacuum', ...
          ott.beam.medium.Medium.DefaultVacuum(), ...
          @(x) isa(x, 'ott.beam.medium.Vacuum'));
      p.parse(varargin{:});

      frequency = 2*pi*p.Results.vacuum.speed ...
          ./p.Results.material.index ./ p.Results.wavelength;

      med = ott.beam.medium.Medium(...
        'material', p.Results.material, ...
        'frequency', frequency, ...
        'vacuum', p.Results.vacuum);
    end
  end

  methods
    function med = Medium(varargin)
      % Construct a new optical medium
      %
      % Usage
      %   medium = Medium(material, frequency, vacuum)
      %
      % Parameters
      %   - material (Material) -- Material properties of medium.
      %
      %   - frequency (numeric) -- Optical (angular) frequency.
      %     Default: ``Medium.DefaultFrequency()``.
      %
      %   - vacuum (Vacuum) -- Vacuum (defines units).
      %     Default: ``Medium.DefaultVacuum()``.

      p = inputParser;
      p.addOptional('material', [], @(x) isa(x, 'ott.beam.medium.Material'));
      p.addOptional('frequency', ...
          ott.beam.medium.Medium.DefaultFrequency(), @isnumeric);
      p.addOptional('vacuum', ...
          ott.beam.medium.Medium.DefaultVacuum(), ...
          @(x) isa(x, 'ott.beam.medium.Vacuum'));
      p.parse(varargin{:});

      med.material = p.Results.material;
      med.frequency = p.Results.frequency;
      med.vacuum = p.Results.vacuum;
    end

    function mat = ott.beam.medium.Material(med)
      % Extract the material property
      mat = med.material;
    end

    function rmat = rdivide(material1, material2)
      % Construct a relative material from two mediums
      %
      % Usage
      %   rmat = material1 ./ material2

      rmat = ott.beam.medium.Relative(material1, material2);
    end

    function b = eq(a, b)
      % Compare two mediums
      %
      % Usage
      %   b = medium1 == medium2

      b = isequaln(a, b);
    end
  end

  methods % Getters/setters
    function val = get.impedance(mat)
      val = sqrt(mat.permeability ./ mat.permittivity);
    end

    function val = get.index(mat)
      val = mat.material.index;
    end

    function mat = set.material(mat, val)
      assert(isa(val, 'ott.beam.medium.Material'), ...
          'material must be a ott.beam.medium.Material');
      mat.material = val;
    end

    function mat = set.frequency(mat, val)
      assert(isnumeric(val) && isscalar(val), ...
          'frequency must be numeric scalar');
      mat.frequency = val;
    end

    function mat = set.vacuum(mat, val)
      assert(isa(val, 'ott.beam.medium.Vacuum'), ...
          'vacuum must be a ott.beam.medium.Vacuum');
      mat.vacuum = val;
    end

    function mat = set.permittivity(mat, val)
      assert(isnumeric(val), 'permittivity must be numeric');
      mat.material.relative_permittivity = val ./ mat.vacuum.permittivity;
    end
    function val = get.permittivity(mat)
      val = mat.material.relative_permittivity .* mat.vacuum.permittivity;
    end

    function mat = set.permeability(mat, val)
      assert(isnumeric(val), 'permeability must be numeric');
      mat.material.relative_permeability = val ./ mat.vacuum.permeability;
    end
    function val = get.permeability(mat)
      val = mat.material.relative_permeability .* mat.vacuum.permeability;
    end

    function val = get.speed(mat)
      val = mat.vacuum.speed ./ mat.material.index;
    end

    function val = get.wavelength(mat)
      val = mat.speed ./ mat.frequency .* (2*pi);
    end
    function val = get.wavelength0(mat)
      val = mat.vacuum.speed ./ mat.frequency .* (2*pi);
    end
  end
end
