classdef InfinitePower
% Defines a constant power property with infinite power.
%
% Properties
%   - power     -- Constant, infinite

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    power
  end

  methods % Getters/setters
    function p = get.power(beam)
      p = Inf;
    end

    function beam = set.power(beam, val)
      error('Beam power is constant and cannot be changed');
    end
  end
end
