classdef Relative < ott.beam.medium.Material
% Describes a ratio between two materials.
% Inherits from :class:`ott.beam.medium.Material`.
%
% Unlike other materials, this material declares properties as ratios
% of two materials.  The relative material can be constructed from two
% materials or mediums using the constructor or using::
%
%   relMaterial = material1 ./ material2
%
% Properties
%   - material1   -- First material
%   - material2   -- Second material
%   - relative_permittivity   -- Relative permittivity between materials
%   - relative_permeability   -- Relative permeability between materials
%   - index                   -- Relative refractive index between materials
%
% Methods
%   - times       -- Dissolve the material or create a new relative material

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    material1
    material2
  end

  properties (Dependent)
    relative_permittivity
    relative_permeability
  end

  methods
    function mat = Relative(material1, material2)
      % Construct a new relative permittivity material
      %
      % Usage
      %   mat = material1 ./ material2
      %
      %   mat = Relative(material1, material2)
      %
      % Parameters
      %   - material1, material2 (Material|other) -- Two materials forming
      %     the relative relationship.  If the objects aren't
      %     :class:`ott.beam.medium.Material` types, attempts to cast
      %     to a Material.

      % Cast materials if needed
      if ~isa(material1, 'ott.beam.medium.Material')
        material1 = ott.beam.medium.Material(material1);
      end
      if ~isa(material2, 'ott.beam.medium.Material')
        material2 = ott.beam.medium.Material(material2);
      end

      % Store properties
      mat.material1 = material1;
      mat.material2 = material2;
    end

    function mat = times(a, b)
      % Dissolve the material or construct a new relative material
      %
      % Usage
      %   mat = material1 .* material2

      if isa(a, 'ott.beam.medium.Relative') && a.material2 == b
        mat = a.material1;    % dissolve
      elseif isa(b, 'ott.beam.medium.Relative') && b.material2 == a
        mat = b.material1;    % dissolve
      else
        mat = a;
        mat.material2 = mat.material2 ./ other;
      end
    end
  end

  methods % Getters/setters

    function val = get.relative_permittivity(mat)
      val = material1.relative_permittivity ./ material2.relative_permittivity;
    end

    function val = get.relative_permeability(mat)
      val = material1.relative_permeability ./ material2.relative_permeability;
    end

    function mat = set.material1(mat, val)
      assert(isa(val, 'ott.beam.medium.Material'), ...
          'material must be a ott.beam.medium.Material');
      mat.material1 = val;
    end

    function mat = set.material2(mat, val)
      assert(isa(val, 'ott.beam.medium.Material'), ...
          'material must be a ott.beam.medium.Material');
      mat.material2 = val;
    end
  end
end
