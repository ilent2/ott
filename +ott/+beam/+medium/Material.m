classdef (Abstract) Material < ott.beam.medium.Medium
% Describes a material medium
%
% Properties
%   - vacuum            -- The background medium (can be a material)
%
% Abstract properties
%   - relative_index    -- Relative refractive index
%   - relative_permittivity -- Relative permittivity
%   - relative_permeability -- Relative permeability
%
% Dependent properties
%   - permittivity      -- Permittivity (units depend on vacuum)
%   - permeability      -- Permeability (units depend on vacuum)
%   - speed             -- Speed (units depend on vacuum)
%   - index             -- Refractive index
%
%   - permittivity0     -- Vacuum medium permittivity
%   - permeability0     -- Vacuum medium permeability
%   - speed0            -- Vacuum medium wave speed
%   - index0            -- Vacuum medium refractive index

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    vacuum
  end

  properties (Abstract)
    relative_index
    relative_permittivity
    relative_permeability
  end

  properties (Dependent)
    permittivity0
    permeability0
    speed0
    index0

    permittivity
    permeability
    speed
    index
  end

  methods
    function mat = Material(vacuum)
      % Abstract constructor for material medium
      %
      % Usage
      %   mat = mat@ott.beam.medium.Material(vacuum)
      %
      % Parameters
      %   - vacuum (ott.beam.medium.Medium)

      mat.vacuum = vacuum;
    end
  end

  methods % Getters/setters
    function mat = set.vacuum(mat, val)
      assert(isa(val, 'ott.beam.medium.Medium'), ...
          'Vacuum must be a ott.beam.medium.Medium');
      mat.vacuum = val;
    end

    function val = get.speed(mat)
      val = mat.speed0 ./ mat.index;
    end
    function val = set.speed(mat, val)
      mat.index = mat.speed0 ./ val;
    end
    function val = get.permittivity(mat)
      val = mat.relative_permittivity * mat.permittivity0;
    end
    function mat = set.permittivity(mat, val)
      mat.relative_permittivity = val ./ mat.permittivity0;
    end
    function val = get.permeability(mat)
      val = mat.relative_permeability * mat.permeability0;
    end
    function mat = set.permeability(mat, val)
      mat.relative_permeability = val ./ mat.permeability0;
    end
    function val = get.index(mat)
      val = mat.relative_index * mat.index0;
    end
    function mat = set.index(mat, val)
      mat.relative_index = val ./ mat.index0;
    end

    function val = get.speed0(mat)
      val = mat.vacuum.speed;
    end
    function val = get.index0(mat)
      val = mat.vacuum.index;
    end
    function val = get.permittivity0(mat)
      val = mat.vacuum.permittivity;
    end
    function val = get.permeability0(mat)
      val = mat.vacuum.permeability;
    end
  end
end
