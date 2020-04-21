classdef MatlabOde < ott.dynamics.Method
% Method wrapping Matlab's ODE solvers
% Inherits from :class:`Method`.
%
% This class provides an interface for the Matlab ode solvers.
% For a list of Matlab ode solvers and when to use them, see
% https://au.mathworks.com/help/matlab/math/choose-an-ode-solver.html
%
% Properties
%   - solver -- Function handle for the Matlab solver

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    solver        % Function handle for the Matlab solver
    options       % Solver options (created with odeset)
    odefun        % Equation describing the problem
    system        % Information about the system
  end

  methods (Hidden)

    function ydash = odeFunStokes(mtd, time, y)
      % Calculate derivatives for stokes step

      % Get position and rotation from input
      position = y(1:3);
      rotation = ott.utils.quaternion2rmatrix(y(4:7));

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

      % Apply drag to force-torque
      dx = inv(drag) * ft;

      % Calculate rotation matrix
      rmatrix = ott.utils.rotation_matrix(dx(4:6));

      % Package output
      ydash = zeros(7, 1);
      ydash(1:3) = dx(1:3);
      ydash(4:7) = ott.utils.rmatrix2quaternion(rmatrix);

    end

    function ydash = odeFunNewtonNoBm(mtd, time, y)
      % Calculate derivatives for inertial step

      ft = mtd.system.force_method(position, rotation, time, drag);

      error('Not yet implemented');

      % TODO: Finish this, how do we do rotations
      ydash = zeros(14, 1);
      ydash(1:7) = y(8:14);
      ydash(8:10) = -drag * y(1:3) ./ mass + ft(1:3) ./ mass;
    end

    function [tout, posOut, rotOut, forceOut] = evaluateInternal(mtd, ...
        tspan, position, rotation)
      % Internal evaluation function, calls Matlab's ode function

      % Convert rotation matrices to quaternions
      quat = ott.utils.rmatrix2quaternion(rotation);

      % Reshape inputs
      y0 = [position(:); quat(:)];

      if isa(mtd.system, 'ott.dynamics.Newtonian')
        y0 = [y0; zeros(size(y0))];
      end

      % Call solver
      [tout, yout] = mtd.solver(mtd.odefun, tspan, y0, mtd.options);

      % Reshape outputs
      posOut = yout(:, 1:3).';
      qOut = yout(:, 4:7).';

      % Convert from quaternions to rotation matrices
      rotOut = ott.utils.quaternion2rmatrix(qOut);

      % Matlab's methods don't seem to support outputting odefun
      % TODO: We could monitor calls to odefun and store outputs
      %   So we don't store too many outputs, we might be able to
      %   trigger storage using odeset.OutputFcn
      forceOut = nan;
    end
  end

  methods
    function mtd = MatlabOde(dynamicsSystem, varargin)
      % Construct a new dynamics simulation method using Matlab's solvers
      %
      % Usage
      %   mtd = MatlabOde(dynamicsSystem, ...)
      %
      % Parameters
      %   - dynamicsSystem -- a :class:`DynamicsSystem` instance.
      %
      % Optional named arguments
      %   - solver (function_handle) -- The name of a Matlab ode solver
      %     or a similar function. For help choosing a solver see
      %     https://au.mathworks.com/help/matlab/math/choose-an-ode-solver.html
      %     Default: ``@ode45``.
      %
      %   - output_function (handle) -- An output function to call after
      %     each simulation step.  For example @odeplot.
      %     Default: ``[]`` (none).

      p = inputParser;
      p.addParameter('solver', @ode45);
      p.addParameter('output_function', []);
      p.parse(varargin{:});

      mtd.system = dynamicsSystem;
      mtd.solver = p.Results.solver;

      % Check for Brownian motion
      if dynamicsSystem.brownian_motion
        warning('ott:dynamics:MatlabOde:no_bm', ...
            'MatlabOde doesn''t support Brownian motion, ignoring');
      end

      if isa(dynamicsSystem, 'ott.dynamics.Stokes')

        % TODO: Can we use Mass with quaternions?
        %mtd.options = odeset('Mass', dynamicsSystem.drag);
        mtd.options = odeset('OutputFcn', p.Results.output_function);

        mtd.odefun = @(t, y) mtd.odeFunStokes(t, y);

      elseif isa(dynamicsSystem, 'ott.dynamics.Newtonian')

        % TODO: Can we use Mass with quaternions?
        %mtd.options = odeset('Mass', dynamicsSystem.mass);
        mtd.options = odeset('OutputFcn', p.Results.output_function);

        mtd.odefun = @(t, y) mtd.odeFunNewton(t, y);

      else
        error('Unknown dynamics system type');
      end
    end
  end

  methods
    function mtd = set.solver(mtd, val)
      % Check solver type
      assert(isa(val, 'function_handle'), ...
          'solver must be a function handle');
      mtd.solver = val;
    end
  end
end
