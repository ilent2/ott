classdef Arbitrary < ott.beam.medium.Material
% Arbitrary magnetic/electronic material
%
% This class declares a material with arbitrary permittivity and
% permeability properties.  These can be isotropic and/or conductive.
%
% Properties
%   - relative_permittivity   -- Scalar or 3x3 matrix
%   - relative_permeability   -- Scalar or 3x3 matrix
%
% Dependent properties
%   - index                 -- Refractive index
%   - isIsotropic           -- True if the material is isotropic
%   - isConductive          -- True if the material is conductive

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    relative_permittivity
    relative_permeability
  end

  methods
    function mat = Arbitrary(varargin)
      % Construct a new arbitrary material definition
      %
      % Usage
      %   mat = Arbitrary(permittivity, permeability)
      %
      % Parameters
      %   - permittivity (numeric) -- Relative permittivity of material.
      %     Can either be a scalar or 3x3 matrix.
      %
      %   - permeability (numeric) -- Relative permeability of material.
      %     Can either be a scalar or 3x3 matrix.
      %
      % Parameters can also be named.

      p = inputParser;
      p.addRequired('permittivity', @isnumeric);
      p.addRequired('permeability', @isnumeric);
      p.parse(varargin{:});

      mat.relative_permittivity = p.Results.permittivity;
      mat.relative_permeability = p.Results.permeability;
    end
  end

  methods % Getters/setters
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

    function mat = set.relative_permeability(mat, val)
      assert(isnumeric(val) && (isscalar(val) ...
          || (ismatrix(val) && all(size(val) == [3, 3]))), ...
          'relative_permeability must be numeric scalar or 3x3 matrix');

      % Convert to scalar if isotropic
      if isdiag(val) && all(diag(val)) == val(1, 1)
        val = val(1, 1);
      end

      mat.relative_permeability = val;
    end
  end
end

