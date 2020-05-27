classdef ZeroPower
% Defines a constant power property with zero power.
%
% Properties
%   - power     -- Constant, zero

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    power
  end

  methods % Getters/setters
    function p = get.power(beam)
      p = 0;
    end

    function beam = set.power(beam, val)
      error('Beam power is constant and cannot be changed');
    end
  end
end
