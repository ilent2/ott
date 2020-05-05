classdef Shape < ott.shapes.Shape
% Describes scattering by a Shape instance.
% Inherits from :class:`ott.shapes.Shape`.
%
% This class describes scattering by Homogeneous shapes with a single
% refractive index value.
%
% Properties
%   - index_relative      -- Refractive index of particle relative to medium
%
% Methods
%   - scatter             -- Calculate scattered rays

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    index_relative
  end

  methods
    function [rray, tray] = scatter(shape, iray)
      % Calculate scattering from an incident ray
      %
      % Usage
      %   [rray, tray] = shape.scatter(iray)

      % TODO

    end
  end
end
