classdef VariablePower
% Adds a variable power property to a Beam.
%
% Properties
%   - power       -- Power property

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    power
  end

  methods % Getters/setters
    function beam = set.power(beam, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0.0, ...
          'Power must be numeric scalar greater than or equal to zero');
      beam.power = val;
    end
  end
end
