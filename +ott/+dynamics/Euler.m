classdef Euler < ott.dynamics.Method
% Simple implementation of Euler's method.
% Inherits from :class:`Method`.
%
% TODO: Early stopping criteria
% TODO: Generate visualisation similar to @odeplot
% TODO: Translation constraints (only allow motion in certain directions)
%
% See also Euler, MatlabOde, MaxStep.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    timestep          % Time-step for the method
    system            % Information about the system
  end

  methods (Hidden)

    function [dx, ft] = step(mtd, position, rotation, time)
      % Calculate a simulation time step
      %
      % Assumes Brownian diffusion

      % Get drag
      if isa(mtd.system.drag, 'function_handle')
        drag = mtd.system.drag(position, rotation, time);
      else
        drag = mtd.system.drag;
      end

      % Apply rotation to drag
      % TODO: Where should this be, in drag calc function?
      drag = drag.rotate(rotation);

      % Calculate force and torque
      ft = mtd.system.force_method(position, rotation, time, drag);

      % Calculate step
      dx = inv(drag) * ft;

      % Add Brownian motion contribution
      if mtd.system.brownian_motion

        % Get temperature
        if isa(mtd.system.temperature, 'function_handle')
          temperature = mtd.system.temperature(position, rotation, time);
        else
          temperature = mtd.system.temperature;
        end

        % Add Brownian correction
        dx = dx + sqrt(2*mtd.system.kb*temperature ...
            ./mtd.timestep)*inv(drag)^(0.5)*randn(size(ft));
      end
    end

    function [tout, posOut, rotOut, forceOut] = evaluateInternal(mtd, ...
        tspan, position, rotation)
      % Simple Euler's method integration

      % Calculate output time steps (possibly skipping last)
      tout = tspan(1):mtd.timestep:tspan(2);

      posOut = zeros(3, length(tout));
      rotOut = zeros(3, 3*length(tout));
      forceOut = zeros(6, length(tout));

      % Store initial position and rotation
      posOut(:, 1) = position;
      rotOut(:, 1:3) = rotation;

      for ii = 2:length(tout)

        % Calculate step
        [dx, ft] = mtd.step(position, rotation, tout(ii));

        % Update position
        position = position + dx(1:3) * mtd.timestep;

        % Update rotation
        rotation = ott.utils.rotation_matrix(dx(4:6)*mtd.timestep) * rotation;

        % Store values
        posOut(:, ii) = position;
        rotOut(:, (1:3) + (ii-1)*3) = rotation;
        forceOut(:, ii-1) = ft;

      end

      % Void final force/torque
      forceOut(:, end) = nan;
    end
  end

  methods
    function mtd = Euler(dynamicsSystem, varargin)
      % Construct a new Eulers method solver instance
      %
      % Usage
      %   mtd = Euler(dynamicsSystem, ...)
      %
      % Optional named arguments
      %   - timestep (numeric) -- Time step size.
      %     Default: 1.0e-3;

      p = inputParser;
      p.addParameter('timestep', 1e-3);
      p.parse(varargin{:});

      mtd.system = dynamicsSystem;
      mtd.timestep = p.Results.timestep;

    end
  end
end

