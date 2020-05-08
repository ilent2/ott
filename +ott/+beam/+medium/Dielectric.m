classdef Dielectric < ott.beam.medium.Material
% Medium describing simple dielectric materials.
%
% In a `Dielectric`, the relative permeability is constant, equal to 1,
% and the refractive index is directly related to the relative permittivity
% by :math:`\sqrt{\epsilon_r}`.
%
% Properties
%   - vacuum                -- Material vacuum
%   - relative_index        -- Refractive index (dependent)
%   - relative_permittivity -- Relative permittivity of medium
%   - relative_permeability -- Relative permeability (constant, 1.0)
%
% Static methods
%   - Water                 -- Material with refractive index of 1.33
%   - Polystyrene           -- Material with refractive index of 1.59
%
% For other properties, see :class:`ott.beam.medium.Material`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    relative_permittivity
  end

  properties (Dependent)
    relative_index
    relative_permeability
  end

  methods (Static)
    function mat = Water(vacuum)
      % Non-conductive dielectric material similar to water at 20 C
      %
      % Refractive index is 1.33.
      %
      % Usage
      %   mat = ott.beam.medium.Dielectric.Water
      %
      %   mat = ott.beam.medium.Dielectric.Water(vacuum)
      %   Specify the vacuum material (default: Unitary)

      if nargin < 1
        vacuum = ott.beam.medium.Vacuum.Unitary;
      end

      mat = ott.beam.medium.Dielectric('vacuum', vacuum, 'index', 1.33);
    end

    function mat = Polystyrene(vacuum)
      % Non-conductive dielectric material similar to polystyrene
      %
      % Refractive index is 1.59.
      %
      % Usage
      %   mat = ott.beam.medium.Dielectric.Polystyrene
      %
      %   mat = ott.beam.medium.Dielectric.Polystyrene(vacuum)
      %   Specify the vacuum material (default: Unitary)

      if nargin < 1
        vacuum = ott.beam.medium.Vacuum.Unitary;
      end

      mat = ott.beam.medium.Dielectric('vacuum', vacuum, 'index', 1.59);
    end
  end

  methods
    function mat = Dielectric(varargin)
      % Construct a new dielectric material
      %
      % Usage
      %   mat = Dielectric(...)
      %
      % Named parameters
      %   - permittivity|relative_permittivity (numeric)
      %   - index|relative_index (numeric)
      %   - speed (numeric)
      %
      %   - vacuum (ott.beam.medium.Medium) -- Medium to use for vacuum.
      %     Default: ``ott.beam.medium.Vacuum.Unitary``
      %
      % Must provide one of permittivity/index/speed.

      p = inputParser;
      p.addParameter('vacuum', ott.beam.medium.Vacuum.Unitary);
      p.addParameter('index', []);
      p.addParameter('relative_index', []);
      p.addParameter('permittivity', []);
      p.addParameter('relative_permittivity', []);
      p.addParameter('speed', []);
      p.parse(varargin{:});

      % Assign vacuum
      mat = mat@ott.beam.medium.Material(p.Results.vacuum);

      % Check refractive index is not ambiguous or undefined
      num_args = isempty(p.Results.index) ...
          + isempty(p.Results.relative_index) ...
          + isempty(p.Results.permittivity) ...
          + isempty(p.Results.relative_permittivity) ...
          + isempty(p.Results.speed);
      assert(num_args == 4, ...
          'ott:beam:medium:Dielectric:wrong_arg_count', ...
          'Must provide one of permittivity/index/speed');

      % Assign refractive index
      if ~isempty(p.Results.index)
        mat.index = p.Results.index;
      elseif ~isempty(p.Results.index_relative)
        mat.index_relative = p.Results.index_relative;
      elseif ~isempty(p.Results.permittivity)
        mat.permittivity = p.Results.permittivity;
      elseif ~isempty(p.Results.permittivity_relative)
        mat.permittivity_relative = p.Results.permittivity_relative;
      elseif ~isempty(p.Results.speed)
        mat.speed = p.Results.speed;
      end
    end
  end

  methods % Getters/setters
    function index = get.relative_index(mat)
      index = sqrt(mat.relative_permittivity);
    end
    function mat = set.relative_index(mat, val)
      assert(isnumeric(val) && isscalar(val), ...
          'relative index must be numeric scalar');
      mat.relative_permittivity = val.^2;
    end

    function val = get.relative_permeability(mat)
      val = 1.0;
    end

    function mat = set.relative_permittivity(mat, val)
      assert(isnumeric(val) && isscalar(val), ...
          'relative permittivity must be numeric scalar');
      mat.relative_permittivity = val;
    end
  end
end

