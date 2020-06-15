classdef InfVolume < ott.shapes.mixin.NoSurfPoints
% Properties associated with infinite volume
% Inherits from :class:`NoSurfPoints`.
%
% Properties
%   - maxRadius
%   - volume

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    maxRadius
    volume
  end

  methods % Getters/setters
    function r = get.maxRadius(shape)
      r = Inf;
    end

    function v = get.volume(shape)
      v = Inf;
    end
  end
end
