classdef NoSurfPoints
% Declares surfPoints method that raise an error.
%
% This class is intended for shapes with infinite surfaces area.
% In these cases it is probably better to choose an iterative method
% or mesh of points based on the specific integration being solved.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function surfPoints(varargin)
      error('No surfPoints method for this shape');
    end
  end
end
