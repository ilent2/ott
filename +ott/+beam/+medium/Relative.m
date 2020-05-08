classdef Relative < ott.beam.medium.Material
% Relative material medium
% Inherits from :class:`ott.beam.medium.Material`.
%
% Properties
%   - material1   -- First material
%   - material2   -- Second material
%
% Methods
%   - times       -- Dissolve the material or create a new relative material
%
% For other properties, see :class:`Material`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    material1
    material2
  end

  properties (Dependent)
    relative_index
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
      %   - material1, material2 (Medium) -- Two mediums forming the
      %     relative permittivity relationship.
      %     If mediums are :class:`Material`, the vacuums should match.
      %
      % Optional named arguments
      %   - vacuum (Medium) -- The material to use for the vacuum.
      %     Default: ``material1.vacuum`` if material is a Material.
      %     Otherwise ``ott.beam.medium.Vacuum.Unitary``.

      p = inputParser;
      p.addParameter('vacuum', ott.beam.medium.Vacuum.Unitary);
      p.parse(varargin{:});

      % Check vacuums match
      if isa(material2, 'ott.beam.medium.Material') ...
          && isa(material1, 'ott.beam.medium.Material')
        assert(material1.vacuum == material2.vacuum, ...
            'vacuum must match when both arguments are Materials');
      end

      % Get vacuum
      if isa(material1, 'ott.beam.medium.Material')
        vacuum = material1.vacuum;
      elseif isa(material2, 'ott.beam.medium.Material')
        vacuum = material2.vacuum;
      else
        vacuum = p.Results.vacuum;
      end

      % TODO: Conversion between vacuums (should be possible?)

      % Construct base
      mat = mat@ott.beam.medium.Material(vacuum);

      % Store materials
      mat.material1 = material1;
      mat.material2 = material2;
    end

    function mat = times(mat, other)
      % Dissolve the material or construct a new relative material
      %
      % Usage
      %   mat = material1 .* material2
      %
      % If `material1.material2 == material2`, returns material1.
      % Otherwise sets `material1.material2` to a new relative
      % material created from `material1.material2` and `material2`.

      if mat.material2 == other
        % Nothing to do
      else
        mat.material2 = mat.material2 ./ other;
      end
    end
  end

  methods % Getters/setters

    function val = get.relative_index(mat)
      val = material1.index ./ material2.index;
    end

    function val = get.relative_permittivity(mat)
      val = material1.permittivity ./ material2.permittivity;
    end

    function val = get.relative_permeability(mat)
      val = material1.permeability ./ material2.permeability;
    end

    function mat = set.material1(mat, val)
      assert(isa(val, 'ott.beam.medium.Medium'), ...
          'material must be a ott.beam.medium.Medium');
      if isa(val, 'ott.beam.medium.Material')
        assert(val.vacuum == mat.vacuum, ...
          'material vacuum must match');
      end
      mat.material1 = val;
    end

    function mat = set.material2(mat, val)
      assert(isa(val, 'ott.beam.medium.Medium'), ...
          'material must be a ott.beam.medium.Medium');
      if isa(val, 'ott.beam.medium.Material')
        assert(val.vacuum == mat.vacuum, ...
          'material vacuum must match');
      end
      mat.material2 = val;
    end
  end
end
