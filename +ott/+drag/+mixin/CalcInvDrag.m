classdef CalcInvDrag
% Calculate inverse drag tensor from forward drag tensor
%
% Abstract properties
%   - forwardInternal
%
% Properties
%   - inverseInternal

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Abstract)
    forwardInternal
  end

  properties (Dependent)
    inverseInternal
  end

  methods % Getters/setters
    function iD = get.inverseInternal(shape)
      iD = inv(shape.forwardInternal);
    end
  end
end

