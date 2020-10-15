classdef VarXySymmetry
% Add property for xySymmetry with value checking
%
% Properties
%   - xySymmetry    -- Variable, default false

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    xySymmetry = false
  end

  methods
    function shape = set.xySymmetry(shape, val)
      assert(islogical(val), 'xySymmetry must be logical');
      shape.xySymmetry = val;
    end
  end
end
