classdef Dielectric < ott.beam.medium.Material
% Medium describing simple dielectric materials.
%
% In a `Dielectric`, the relative permeability is constant, equal to 1,
% and the refractive index is directly related to the relative permittivity
% by :math:`n = \sqrt{\epsilon_r}`.
%
% Properties
%   - relative_permittivity -- Relative permittivity of medium
%   - relative_permeability -- Relative permeability (constant, 1.0)
%
% Dependent properties
%   - index                 -- Refractive index
%   - isIsotropic           -- True if the material is isotropic
%   - isConductive          -- True if for electrically conductive materials
%
% Static methods
%   - FromIndex             -- Construct medium from refractive index
%
% See also :class:`ott.beam.medium.Generic` for example materials.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    relative_permittivity
  end

  properties (Dependent)
    relative_permeability
  end

  methods (Static)
    function mat = FromIndex(index)
      % Construct material from refractive index
      %
      % Usage
      %   mat = Dielectric.FromIndex(index)

      assert(isnumeric(index), 'index must be numeric');
      mat = ott.beam.medium.Dielectric(index.^2);
    end
  end

  methods
    function mat = Dielectric(varargin)
      % Construct a new dielectric material
      %
      % Usage
      %   mat = Dielectric(permittivity)
      %
      % Parameters
      %   - permittivity (numeric) -- Relative permittivity of material.
      %     Can either be a scalar or 3x3 matrix.
      %
      % See also :meth:`FromIndex` for a constructor from refractive index.

      p = inputParser;
      p.addRequired('permittivity', @isnumeric);
      p.parse(varargin{:});

      mat.relative_permittivity = p.Results.permittivity;
    end
  end

  methods % Getters/setters
    function val = get.relative_permeability(mat)
      val = 1.0;
    end

    function mat = set.relative_permittivity(mat, val)
      assert(isnumeric(val) && (isscalar(val) ...
          || (ismatrix(val) && all(size(val) == [3, 3]))), ...
          'relative_permittivity must be numeric scalar or 3x3 matrix');

      % Convert to scalar if isotropic
      if isdiag(val) && all(diag(val)) == val(1, 1)
        val = val(1, 1);
      end

      mat.relative_permittivity = val;
    end
  end
end

