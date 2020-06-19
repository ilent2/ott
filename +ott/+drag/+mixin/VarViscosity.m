classdef VarViscosity
% Declares a variable viscosity field.
%
% No constructor is implemented.  No default value is set.
%
% Properties
%   - viscosity

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    viscosity
  end

  methods % Getters/setters
    function drag = set.viscosity(drag, val)
      assert(isnumeric(val) && isscalar(val), ...
        'viscosity must be numeric scalar');
      drag.viscosity = val;
    end
  end
end

