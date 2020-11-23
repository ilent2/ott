classdef WallEffect < ott.dynamics.Dynamics
% Dynamics simulation with a hydrodynamic wall effect.
% Inherits from :class:`Dynamics`.
%
% This class implements a dynamics simulation where the drag tensor is
% re-calculated at each step of the simulation using a wall effect method.
%
% This class may change in future to include an option for optical
% interaction, for now it only  includes hydrodynamic interaction.
%
% Properties
%   - dragMethod   -- Method to use to calculate drag.
%
% See base :class:`Dynamics` for additional properties/methods.

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    dragMethod      % Method used to calculate drag (function_handle)
  end

  methods
    function sim = WallEffect(varargin)
      % Create a wall effect dynamics instance
      %
      % Usage
      %   sim = WallEffect(...)
      %
      % Parameters
      %   - dragMethod (function_handle) -- A function handle specifying
      %     a drag calculation method.  Default::
      %
      %       @(shape) ott.drag.StokesSphereWall.FromShape(
      %           shape, 'viscosity', 0.001)``.
      %
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
      p.addParameter('dragMethod', ...
        @(shape) ott.drag.StokesSphereWall.FromShape(shape, ...
        'viscosity', 0.001));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      sim = sim@ott.dynamics.Dynamics(unmatched{:});
      sim.dragMethod = p.Results.dragMethod;
    end

    function i0 = setupSimulation(sim, ~, x, R)
      % Called to setup the simulation
      %
      % This is exactly the same as Isolated.setupSimulation

      assert(~isempty(sim.particle) && ~isempty(sim.particle.shape), ...
          'A shape instance is required for this method');
      shape = sim.particle.shape;
      assert(numel(shape) == 2 && (isa(shape(1), 'ott.shape.Plane') ...
          || isa(shape(2), 'ott.shape.Plane')), ...
          'shape must be a 2 element array containing a plane');

      % Flip shape array if required
      if isa(shape(1), 'ott.shape.Plane')
        sim.particle.shape = flip(sim.particle.shape);
      end

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
      % This is slightly modified from Isolated.simulationStep

      persistent internalStep optforce odx

      if t < 0
        % Setup the simulation

        kb = sim.boltzmannConstant;
        T = sim.temperature;
        dt = sim.timeStep;
        bmterm = sqrt(2*kb*T./dt);

        % Get the default force method
        optforce = sim.opticalForceMethod;

        % Two different loops depending on if we need inertia
        if (isempty(sim.particle.mass) && isempty(sim.particle.moment)) ...
            || (sim.particle.mass == 0 && sim.particle.moment == 0)

          % No mass, no inertia

          internalStep = @(gft, drag) inv(drag)*gft ...
              + bmterm*drag^(-1/2)*randn(6, 1); %#ok<MINV>

        else

          % Mass provided, include inertia

          % Get mass
          mass = sim.particle.mass;
          moment = sim.particle.moment;
          if isempty(mass), mass = 0.0; end
          if isempty(moment), moment = 0.0; end
          massmoment = [[1;1;1]*mass; [1;1;1]*moment];
          diagMoment = diag(massmoment);

          internalStep = @(gft, drag) inv(diagMoment + dt*drag) ...
              * (odx.*massmoment + gft*dt ...
              + dt*bmterm*drag^(1/2)*randn(6, 1)); %#ok<MINV>

          % Initial derivatives
          % TODO: Add initial rotation derivative
          odx = zeros(6, 1);
          odx(1:3) = (x(:, 2) - x(:, 1))./dt;

        end
        
        return;

      end
      
      % Update particle position
      % We cant make shape persistent since we use it for the figure
      sim.particle.shape(1).position = x;
      
      % Check if particle has collided with wall
      if sim.hasCollided
        dx = nan(3, 1);
        return;
      end

      % Calculate force/torque
      gft = -optforce(t, x, R);

      % Calculate drag at this position
      drag = sim.dragMethod(sim.particle.shape);

      % Calculate step
      dx = internalStep(gft, drag);
      odx = dx;

    end
    
    function b = hasCollided(sim)
      % Check if the particle has collided with the wall
      %
      % Usage
      %   b = sim.hasCollided()
      
      norm = sim.particle.shape(2).normal;
      dist = dot(norm, sim.particle.shape(1).position ...
        - sim.particle.shape(2).position);
      b = dist <= sim.particle.shape(1).radius;
      
    end
  end

  methods % Getters/setters
    function sim = set.dragMethod(sim, val)
      assert(isa(val, 'function_handle') && nargin(val) == 1, ...
          'dragMethod must be a function handle taking 1 argument');
      sim.dragMethod = val;
    end
  end
end

