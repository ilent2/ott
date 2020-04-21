classdef (Abstract) Method
% Base class for optical tweezers dynamics methods.
%
% Static methods
%   - simple -- Attempt to automatically choose a good method.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function mtd = simple(dynamicsSystem)
      % Attempts to automatically choose a method based on the system
      %
      % Usage
      %   mtd = Method.simple(dynamicsSystem)
      %
      % Parameters
      %   - dynamicsSystem -- a :class:`DynamicsSystem` instance.
      
      % TODO: Smart decisions about methods
      mtd = ott.dynamics.Euler(dynamicsSystem);
    end
  end

  methods (Abstract, Hidden)
    % Method which must be implemented by derived classes
    % This allows us to enforce the interface for the evaluate function
    evaluateInternal
  end

  methods
    function [tout, posOut, rotOut, forceOut] = evaluate(mtd, ...
        tspan, position, rotation)
      % Evaluate the initial value problem for the specified time period
      %
      % Usage
      %   [tout, posOut, rotOut, forceOut] = mtd.evaluate(tspan, ...
      %       position, rotation)
      %   Simulates the dynamics and outputs the recorded steps, positions,
      %   rotations and calculated forces.
      %
      % Although forces aren't a part of the ode solution, they are
      % typically expensive to calculate so we store them too.
      %
      % Parameters
      %   - tspan (2xnumeric) -- Time span for simulation
      %   - position (numeric) -- Initial particle position
      %   - rotation (numeric) -- Initial particle orientation

      [tout, posOut, rotOut, forceOut] = mtd.evaluateInternal(...
          tspan, position, rotation);
    end
  end
end
