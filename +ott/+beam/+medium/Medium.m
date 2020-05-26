classdef (Abstract) Medium
% Abstract base class for describing optical media
%
% Abstract properties
%   - permittivity  -- Permittivity of medium
%   - permeability  -- Permeability of medium
%   - speed         -- Wave speed in medium
%   - index         -- Refractive index of medium
%
% Dependent properties
%   - impedance     -- Impedance of medium
%
% Methods
%   - rdivide       -- Create a relative material
%
% See also :class:`ott.beam.medium.Vacuum`, :class:`ott.beam.medium.Material`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties (Abstract)
    permittivity
    permeability
    speed
    index
  end

  properties (Dependent)
    impedance
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list

      % Nothing to do
    end
  end

  methods
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
      % TODO: Should this really be a dependent property
      % TODO: Should this be matrix multiplication/sqrt?
      val = sqrt(mat.permeability ./ mat.permittivity);
    end
  end
end
