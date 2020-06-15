classdef VarStarShaped
% Add property for starShaped with value checking
%
% Properties
%   - starShaped    -- Variable, default false

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    starShaped = false
  end

  methods
    function shape = set.starShaped(shape, val)
      assert(islogical(val) && isscalar(val), ...
          'starShaped must be scalar logical');
      shape.starShaped = val;
    end
  end
end
