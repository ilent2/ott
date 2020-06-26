classdef (Abstract) Material
% Description of an optical material (no units or frequency)
%
% Properties (Abstract)
%   - relative_permittivity    -- Relative permittivity of medium
%   - relative_permeability    -- Relative permeability of medium
%
% Dependent properties
%   - index           -- Refractive index
%   - isIsotropic     -- True if the material is isotropic
%   - isConductive    -- True if the material is conductive
%   - isElectric      -- True if the permittivity is not 1
%   - isMgnetic       -- True if the permeability is not 1
%
% Methods
%   - rdivide       -- Create a relative material

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties (Abstract)
    relative_permittivity
    relative_permeability
  end

  properties (Dependent)
    index             % Refractive index (or relative refractive index)
    isIsotropic       % True if isotropic
    isConductive      % True if conductive
    isElectric        % True if the permittivity is not 1
    isMgnetic         % True if the permeability is not 1
  end

  methods
    function rmat = rdivide(material1, material2)
      % Construct a relative material from two mediums
      %
      % Usage
      %   relMaterial = material1 ./ material2

      rmat = ott.beam.medium.Relative(material1, material2);
    end
  end

  methods % Getters/setters
    function val = get.index(mat)
      val = sqrt(mat.relative_permittivity .* mat.relative_permeability);
    end

    function val = get.isIsotropic(mat)
      val = isscalar(mat.relative_permittivity) ...
          & isscalar(mat.relative_permeability);
    end

    function val = get.isConductive(mat)
      val = ~isreal(mat.relative_permittivity);
    end

    function val = get.isElectric(mat)
      val = mat.relative_permittivity ~= 1.0;
    end

    function val = get.isMagnetic(mat)
      val = mat.relative_permeability ~= 1.0;
    end
  end
end
