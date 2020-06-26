classdef VariablePower
% Adds a variable power property to a Beam.
%
% Doesn't provide a constructor, uses properties.Beam.
%
% Properties
%   - power       -- Power property
%
% Static methods
%   - likeProperties    -- Adds power to the properties list

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    power
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      if isa(other, 'ott.beam.prperties.VariablePower')
        args = ott.utils.addDefaultParameter('power', other.power, args);
      end
    end
  end

  methods % Getters/setters
    function beam = set.power(beam, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0.0, ...
          'Power must be numeric scalar greater than or equal to zero');
      beam.power = val;
    end
  end
end
