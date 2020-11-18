classdef Isolated < ott.dynamics.Dynamics
% Setup a dynamics simulation of an isolated particle.
% Inherits from :class:`Dynamics`.
%
% This class supports simulating isolated single particles in (or outside)
% an optical trap.  Drag, T-matrix and beam properties all need to be
% calculated before running the simulation.
%
% For further details/properties/methods, see :class:`Dynamics`.

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function sim = Isolated(varargin)
      % Create a new dynamics instance, storing properties/setting defaults
      %
      % Usage
      %   sim = Isolated(...)
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

      % Pass all to base
      sim = sim@ott.dynamics.Dynamics(varargin{:});
    end

    function i0 = setupSimulation(sim, ~, x, R)
      % Called to setup the simulation

      % Two different loops depending on if we need inertia
      if (isempty(sim.particle.mass) && isempty(sim.particle.moment)) ...
          || (sim.particle.mass == 0 && sim.particle.moment == 0)

        % No mass, no inertia
        i0 = 2;

        % Call step function for setup
        sim.simulationStep(-1, [], []);

      else

        % Mass provided, include inertia
        i0 = 3;

        % Call step function for setup
        sim.simulationStep(-1, x(:, 1:2), R(:, 1:6));

      end

    end

    function dx = simulationStep(sim, t, x, R)
      % Called to advance the simulation
      %
      % This function is implemented without using class properties,
      % in the past we have seen significant overhead when using properties.
      % It is hoped by using persistent variables we can avoid this overhead.

      persistent internalStep optforce odx

      if t < 0
        % Setup the simulation

        kb = sim.boltzmannConstant;
        dt = sim.timeStep;

        % Get the default force method
        optforce = sim.opticalForceMethod;

        % Two different loops depending on if we need inertia
        if (isempty(sim.particle.mass) && isempty(sim.particle.moment)) ...
            || (sim.particle.mass == 0 && sim.particle.moment == 0)

          % No mass, no inertia

          if ~isempty(sim.particle.drag)
            invGamma = inv(sim.particle.drag);
            bmterm = sqrt(2*kb*sim.temperature./dt) * invGamma^(1/2);
            internalStep = @(gft) invGamma*gft + bmterm*randn(6, 1); %#ok<MINV>
          else
            internalStep = @(gft) gft;
          end

        else

          % Mass provided, include inertia

          % Get mass
          mass = sim.particle.mass;
          moment = sim.particle.moment;
          if isempty(mass), mass = 0.0; end
          if isempty(moment), moment = 0.0; end

          if ~isempty(sim.particle.drag)
            drag = sim.particle.drag;
            T = sim.temperature;
          else
            drag = eye(6);
            T = 0;
          end
          massmoment = [[1;1;1]*mass; [1;1;1]*moment];
          denom = inv(diag(massmoment) + dt * drag);
          bmterm = sqrt(2*kb*T./dt) * drag^(1/2);

          internalStep = @(gft) denom * (odx.*massmoment + gft*dt ...
              + dt*bmterm*randn(6, 1)); %#ok<MINV>

          % Initial derivatives
          % TODO: Add initial rotation derivative
          odx = zeros(6, 1);
          odx(1:3) = (x(:, 2) - x(:, 1))./dt;

        end

        return;
      end

      % Calculate force/torque
      gft = -optforce(t, x, R);

      % Calculate step
      dx = internalStep(gft);
      odx = dx;

    end
  end
end
