classdef VarZRotSymmetry
% Add property for zRotSymmetry with value checking
%
% Properties
%   - zRotSymmetry    -- Variable, default 1

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    zRotSymmetry = 1
  end

  methods
    function shape = set.zRotSymmetry(shape, val)
      assert(isscalar(val), 'zRotSymmetry must be scalar');
      shape.zRotSymmetry = val;
    end
  end
end
