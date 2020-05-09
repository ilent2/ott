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
%   - impedance         -- Impedance of the medium
%   - permittivity      -- Permittivity (units depend on vacuum)
%   - permeability      -- Permeability (units depend on vacuum)
%   - speed             -- Speed (units depend on vacuum)
%   - index             -- Refractive index
%
%   - permittivity0     -- Vacuum medium permittivity
%   - permeability0     -- Vacuum medium permeability
%   - speed0            -- Vacuum medium wave speed
%   - index0            -- Vacuum medium refractive index
%
% Static methods
%   - Simple            -- Construct a material from a property list

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

  methods (Static)
    function mat = Simple(varargin)
      % Construct a new material instance from a list of properties
      %
      % Usage
      %   Material.Simple(...)
      %
      % Named parameters
      %   - vacuum (ott.beam.medium.Medium) -- Medium to use as the
      %     vacuum material.  Can be another material.
      %     Default: ``ott.beam.medium.Vacuum.Unitary``.
      %
      %   - like (ott.beam.media.Medium|Material) -- If present, uses
      %     this material or medium for default values.
      %
      % The vacuum is constructed from :meth:`Vacuum.Simple`.
      % All unmatched parameters along with `vacuum` are passed to
      % this function.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('vacuum', []);
      p.addParameter('like', []);
      p.addParameter('permittivity', []);
      p.addParameter('permeability', []);
      p.addParameter('speed', []);
      p.addParameter('index', []);
      p.addParameter('relative_index', []);
      p.addParameter('relative_permittivity', []);
      p.addParameter('relative_permeability', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get default vacuum parameters
      default_vacuum = ott.beam.medium.Vacuum.Unitary;
      if ~isempty(p.Results.like)
        if isa(p.Results.like, 'ott.beam.medium.Material')
          default_vacuum = p.Results.like.vacuum;
        end
      end

      % Get vacuum from inputs or default
      vacuum = p.Results.vacuum;
      if isempty(vacuum)
        vacuum = default_vacuum;
      end

      % Construct vacuum (consuming remaining arguments)
      vacuum = ott.beam.medium.Vacuum.Simple(...
          'like', vacuum, unmatched{:});

      % Get default material parameters
      default_permittivity = 1.0;
      default_permeability = 1.0;
      default_index = 1.0;
      if ~isempty(p.Results.like)
        default_permittivity = p.Results.like.permittivity ...
            ./ vacuum.permittivity;
        default_permeability = p.Results.like.permeability ...
            ./ vacuum.permeability;
        default_index = p.Results.like.index ./ vacuum.index;
      end

      % Check for extra arguments
      assert(isempty(p.Results.relative_permittivity) ...
          || isempty(p.Results.permittivity), ...
          'Too many permittivity arguments supplied.');
      assert(isempty(p.Results.relative_permeability) ...
          || isempty(p.Results.permeability), ...
          'Too many permeability arguments supplied.');
      assert(isempty(p.Results.relative_index) ...
          || isempty(p.Results.index), ...
          'Too many index/speed arguments supplied.');

      % Convert from non-relative inputs to relative inputs
      udef_permittivity = false;
      udef_permeability = false;
      udef_index = false;
      rel_permittivity = p.Results.relative_permittivity;
      rel_permeability = p.Results.relative_permeability;
      rel_index = p.Results.relative_index;
      if isempty(rel_permittivity)
        if ~isempty(p.Results.permittivity)
          rel_permittivity = p.Results.permittivity
          rel_permittivity = rel_permittivity ./ vacuum.permittivity;
        else
          rel_permittivity = default_permittivity;
          udef_permittivity = true;
        end
      end
      if isempty(rel_index)
        if ~isempty(p.Results.index)
          rel_index = p.Results.index;
          rel_index = rel_index ./ vacuum.index;
        else
          rel_index = default_index;
          udef_index = true;
        end
      end
      if isempty(rel_permeability)
        if ~isempty(p.Results.permeability)
          rel_permeability = p.Results.permeability
          rel_permeability = rel_permeability ./ vacuum.permeability;
        else
          rel_permeability = default_permeability;
          udef_permeability = true;
        end
      end

      % Construct material using appropriate constructor
      tol = 1.0e-6;
      if abs(rel_permeability - 1.0) < tol
        if ~udef_index
          mat = ott.beam.medium.Dielectric('vacuum', vacuum, ...
              'relative_index', rel_index);
        else
          mat = ott.beam.medium.Dielectric('vacuum', vacuum, ...
              'relative_permittivity', rel_permittivity);
        end
      else
        % TODO: Add support for magnetic materials
        error('Only dielectric mediums supported for not');
      end
    end
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
