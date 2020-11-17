classdef Dynamics
% Simple dynamics simulations with a fixed time step.
%
% This class stores properties related to dynamics simulations and contains
% methods for running simple dynamics simulations with fixed step sizes.
% Dynamics is described by
%
%     m\ddot{x} = F_{Optical} + F_{Drag} + F_{BM}
%
% If ``temperature`` is zero, the Brownian motion term is omitted. Drag
% is only included if the ``particle.drag`` property is set.  Optical
% force is only included when the ``particle.tmatrix`` property is set.
% Inertia is included when ``particle.mass`` is set.
%
% In a future release, this class will change/move.  It would be nice to
% support other variable step size dynamics simulations.
%
% Properties
%   - beam        -- A optical tweezers toolbox beam
%   - particle    -- A particle which scatters the beam
%   - temperature -- Temperature for Brownian motion term [K]
%   - timeStep    -- Simulation time step [s]
%
% Methods
%   - simulate    -- Run dynamics simulation
%
% Static methods
%   - rotation_matrix -- Convert from torque to rotation matrix

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    beam            % A optical tweezers toolbox beam
    particle        % A particle which scatters the beam
    temperature     % Temperature for Brownian motion term [K]
    timeStep        % Simulation time step
  end

  methods (Static)
    function R = rotation_matrix(rot_axis,rot_angle)
      % Calculates rotation matrix using Euler-Rodrigues formula.
      %
      % Usage
      %   R = ott.tools.Dynamics.rotation_matrix( axis, angle ) --
      %   calculates the rotation about axis by angle (in radians).
      %
      %   R = rotation_matrix( axis_angle ) -- calculates the rotation
      %   about vector axis_angle, the angle is specified as the length
      %   of the vector (in radians).

      if nargin < 2
        rot_angle = vecnorm(rot_axis);
        rot_axis = rot_axis ./ rot_angle;
      end

      assert(isnumeric(rot_angle) && isscalar(rot_angle), ...
          'rot_angle must be numeric scalar');
      assert(isnumeric(rot_axis) && isvector(rot_axis) ...
        && numel(rot_axis) == 3, ...
        'rot_axis should be 3 element numeric vector');

      % Check for small angle
      tol = 1.0e-8;
      if abs(rot_angle) < tol
        R = eye(3);
        return;
      end

      % Check axis is unit vector
      if abs(vecnorm(rot_axis)-1) > 2*eps(1)
        warning('ott:utils:rotation_matrix:unit_vector', ...
          'Rotation axis (with two inputs) should be unit vector, rescaling');
        rot_axis = rot_axis ./ vecnorm(rot_axis);
      end

      % Euler-Rodrigues formula, might be unstable for small angles :(

      ux = [ 0 -rot_axis(3) rot_axis(2);
          rot_axis(3) 0 -rot_axis(1);
          -rot_axis(2) rot_axis(1) 0 ];

      R = eye(3) + sin(rot_angle) * ux + ( 1 - cos(rot_angle) ) * ux * ux;
    end
  end

  methods
    function sim = Dynamics(varargin)
      % Create a new dynamics instance, storing properties/setting defaults
      %
      % Usage
      %   sim = Dynamics(...)
      %
      % Optional named parameters
      %   - temperature (numeric) -- Temperature [K] for Brownian motion
      %     simulations.  Set to 0 for no Brownian motion.  Default: ``300``.
      %
      %   - beam (ott.beam.Beam) -- Beam for optical force calculation.
      %     Default: ``ott.beam.Empty()``
      %
      %   - particle (ott.particle.Particle) -- Description of the
      %     particle.  Default: ``ott.particle.Fixed()``.
      %
      %   - timeStep (numeric) -- Simulation time step.  You may need
      %     to make this smaller, depending on the particle/beam/drag.
      %     Default: ``1.0e-4``.

      p = inputParser;
      p.addParameter('temperature', 300.0);
      p.addParameter('beam', ott.beam.Empty());
      p.addParameter('particle', ott.particle.Fixed());
      p.addParameter('timeStep', 1e-4);
      p.parse(varargin{:});

      sim.beam = p.Results.beam;
      sim.particle = p.Results.particle;
      sim.temperature = p.Results.temperature;
      sim.timeStep = p.Results.timeStep;
    end

    function [t, x, R] = simulate(sim, varargin)
      % Run the dynamics simulation
      %
      % Usage
      %   particle.simulate(...) -- Runs the dynamics simulation, showing
      %   the results in the current figure window (along with a stop
      %   button).
      %
      %   [t, x, R] = particle.simulate(time, ...) -- Runs the simulation
      %   for a fixed duration, returning the times, positions and
      %   rotations of the particle.  No figure is displayed by default.
      %
      % Returns
      %   - t (N numeric) -- Simulation times
      %   - x (3xN numeric) -- Particle positions.
      %   - R (3x3N numeric) -- Particle orientations.
      %
      % Parameters
      %   - time (numeric) -- Amount of time to simulate particle for.
      %     If time isn't divisible by the time step, rounds up.
      %
      % Optional named parameters
      %   - position (3x1 | 3x2 numeric) -- Initial particle position.
      %     Default: ``particle.position``.  For simulations with inertia,
      %     specify the last two particle positions.
      %
      %   - rotation (3x3 | 3x6 numeric) -- Initial particle orientation.
      %     Default: ``particle.rotation``.  For simulations with inertia,
      %     specify the last two particle orientations.
      %
      %   - plot_axes (axes | []) -- Handle to the plot axes or empty
      %     for no plot.  Default: ``gca()`` if no time specified,
      %     otherwise ``[]``.
      %
      %   - outputRate (numeric) -- How many seconds to wait between
      %     updating the figure.  Default: ``1``.

      p = inputParser;
      p.addOptional('time', [], @isnumeric);
      p.addParameter('position', sim.particle.position);
      p.addParameter('rotation', sim.particle.rotation);
      p.addParameter('plot_axes', []);
      p.addParameter('outputRate', 1);
      p.parse(varargin{:});

      time = p.Results.time;
      assert(isempty(time) || (isnumeric(time) ...
          && isscalar(time) && time > 0), ...
          'time must be empty or positive numeric scalar');

      % Setup the figure if required
      our_axes = p.Results.plot_axes;
      if isempty(time) || ~isempty(our_axes)
        plotData = sim.setupAxes(our_axes, ...
            p.Results.outputRate, sim.particle.shape);
      else
        plotData = sim.setupAxes();   % Create empty plotData
      end

      % Calculate number of time steps
      numSteps = ceil(time ./ sim.timeStep);
      t = 0:sim.timeStep:numSteps*sim.timeStep;

      % Allocate memory for output
      x = zeros(3, numel(t));
      R = zeros(3, 3*numel(t));

      % Get initial position/rotation
      if size(p.Results.position, 2) == 2
        x(:, 1:2) = p.Results.position;
      else
        x(:, 1:2) = repmat(p.Results.position, 1, 2);
      end
      if size(p.Results.rotation, 2) == 6
        R(:, 1:6) = p.Results.rotation;
      else
        R(:, 1:6) = repmat(p.Results.rotation, 1, 2);
      end

      kb = 1.3806e-23;    % Boltzmann constant [m2.kg/s2/K]
      dt = sim.timeStep;

      % Optical force calculation function
      if ~isempty(sim.particle.tmatrix)
        optforce = @(x, R) forcetorque(sim.beam, sim.particle, ...
            'position', x, 'rotation', R);
      else
        optforce = @(~, ~) deal(zeros(3, 1), zeros(3, 1));
      end

      % Two different loops depending on if we need inertia
      if (isempty(sim.particle.mass) && isempty(sim.particle.moment)) ...
          || (sim.particle.mass == 0 && sim.particle.moment == 0)

        % No mass, no inertia

        if ~isempty(sim.particle.drag)
          invGamma = inv(sim.particle.drag);
        else
          invGamma = eye(6);
        end
        bmterm = sqrt(2*kb*sim.temperature) * invGamma^(1/2) ...
          ./ sqrt(sim.timeStep);

        for ii = 2:numel(t)

          % Get current position/rotation
          xc = x(:, ii-1);
          Rc = R(:, (1:3) + (ii-2)*3);

          % Calculate force/torque
          [fo, to] = optforce(xc, Rc);

          % Convert to position units and add BM
          ft = -invGamma * [fo; to] + bmterm * randn(6, 1);

          % Update position/rotation
          x(:, ii) = xc + ft(1:3)*dt;
          R(:, (1:3) + (ii-1)*3) = sim.rotation_matrix(ft(4:6)*dt)*Rc;

          % Update plot
          plotData = sim.updatePlot(plotData, ...
              x(:, ii), R(:, (1:3) + (ii-1)*3));

          % Check if we should stop
          if ~plotData.running
            break;
          end
        end

      else

        % With mass, get a local copy
        mass = sim.particle.mass;
        moment = sim.particle.moment;
        if ~isempty(sim.particle.drag)
          drag = sim.particle.drag;
        else
          drag = eye(3);
        end
        massmoment = [[1;1;1]*mass; [1;1;1]*moment];
        denom = inv(diag(massmoment) + dt * drag);
        bmterm = sqrt(2*kb*sim.temperature./dt) * drag^(1/2);
        
        % Initial derivatives
        % TODO: Add initial rotation derivative
        dx = zeros(6, 1);
        dx(1:3) = (x(:, 2) - x(:, 1))./dt;

        for ii = 3:numel(t)

          % Calculate optical force/torque
          [fo, to] = optforce(x(:, ii-1), R(:, (1:3) + (ii-2)*3));
          gft = -[fo; to];
          
          % Calculate new derivative
          dx = denom * (dx.*massmoment + gft.*dt + dt*bmterm*randn(6, 1)); %#ok<MINV>

          % Update position/rotation
          x(:, ii) = x(:, ii-1) + dx(1:3)*dt;
          R(:, (1:3) + (ii-1)*3) = sim.rotation_matrix(dx(4:6)*dt)*R(:, (1:3) + (ii-2)*3);

          % Update plot
          plotData = sim.updatePlot(plotData, ...
              x(:, ii), R(:, (1:3) + (ii-1)*3));

          % Check if we should stop
          if ~plotData.running
            break;
          end
        end

      end

      % Remove extra entries in t/x/R
      if ~plotData.running
        t(ii+1:end) = [];
        x(:, ii+1:end) = [];
        R(:, ii*3+1:end) = [];
      end

      % Only return results if requested
      if nargout == 0
        clear t x R
      end
    end
  end

  methods (Hidden, Static)
    function plotData = setupAxes(oaxes, outputRate, shape)
      % Setup a axes for visualising the simulation

      if nargin == 0
        plotData = struct('running', true, 'axes', []);
        return;
      end

      % setup axes if needed
      if isempty(oaxes)
        oaxes = gca();
      end

      % Add stop button
      stopButton = uicontrol('Parent', oaxes.Parent, ...
          'Style','pushbutton', 'String','Stop', ...
          'Units','normalized', ...
          'Position', [0.0046 0.0071 0.1329 0.0648], ...
          'Visible','on', 'FontSize', 10, ...
          'Callback', @buttonCallback);

      % Get or generate default shape
      if isempty(shape)
        shape = ott.shape.Sphere(1.0e-6);
      end

      % Generate patch data
      opatch = shape.surf('axes', oaxes);

      % Setup plot data
      plotData = struct();
      plotData.axes = oaxes;
      plotData.patch = opatch;
      plotData.patchVertices = opatch.Vertices.';
      plotData.running = true;
      plotData.stopButton = stopButton;
      plotData.time = now;
      plotData.outputRate = outputRate;

      function buttonCallback(~, ~)
        stopButton.UserData = 'stop';
      end
    end

    function plotData = updatePlot(plotData, x, R)
      % Update the plot

      % Check if we need to do anything
      if isempty(plotData.axes)
        return;
      end

      % Check if the user closed the figure
      if ~ishandle(plotData.axes)
        plotData.running = false;
        return;
      end

      % Check if we should update the figure content
      if (now - plotData.time) > plotData.outputRate/86400

        % Update patch position
        plotData.patch.Vertices = (R*plotData.patchVertices + x).';

        % Update content and run callback functions
        drawnow;

        % Check if we have been asked to exit
        if ~isempty(plotData.stopButton.UserData)
          plotData.running = false;
        end

        plotData.time = now;
      end
    end
  end

  methods % Getters/setters
    function sim = set.beam(sim, val)
      assert(isa(val, 'ott.beam.Beam') && isscalar(val), ...
          ['beam must be a single ott.beam.Beam instance', newline, ...
          'Use `ott.beam.Array` for arrays of beams']);
      sim.beam = val;
    end

    function sim = set.particle(sim, val)
      assert(isa(val, 'ott.particle.Particle') && isscalar(val), ...
          'particle must be a single ott.particle.Particle instance');
      sim.particle = val;
    end

    function sim = set.temperature(sim, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0, ...
          'temperature must be positive numeric scalar');
      sim.temperature = val;
    end

    function sim = set.timeStep(sim, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'timeStep must be positive numeric scalar');
      sim.timeStep = val;
    end
  end
end

