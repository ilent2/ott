classdef Euler < ott.dynamics.Method
% Simple implementation of Euler's method.
% Inherits from :class:`Method`.
%
% See also Euler, MatlabOde, MaxStep.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    timestep          % Time-step for the method
  end

  methods
    function step(mtd)

      % What about inertia?
      % What about stochastic terms?

      x = x + v * mtd.teimstep;

    end
  end

end

