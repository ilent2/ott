classdef (Abstract) DynamicsSystem
% Base class for dynamics system
%
% These classes are property collections and identifiers used by the
% other dynamics methods.  The actual DEs are implemented as static
% methods of the dynamics classes.
%
% Static methods
%   - simple -- Automatically attempt to choose appropriate system
%
% Properties
%   - force_method          -- Force calculation method
%   - drag                  -- Drag tensor
%   - temperature           -- Temperature for brownian motion
%   - brownian_motion       -- If brownian motion is used

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    force_method           % Force calculation method
    drag                   % Drag tensor
    temperature            % Temperature for brownian motion
    brownian_motion        % If brownian motion is used
  end

  methods (Static)

    function sys = simple(forceMethod, varargin)
      % Automatically choose a dynamics system based on inputs.
      %
      % If drag and mass are methods, the method selection will
      % use the position, rotation and time optional arguments.
      % This may result in the method selection being sub-optimal.
      %
      % Usage
      %   sys = DynamicsSystem.simple(forceMethod,, ...)
      %   Constructs a new dynamics system for the forceMethod.
      %
      % Optional named arguments
      %   - drag -- Simulation drag or function handle.
      %     Default: ``1.0``.
      %
      %   - mass (numeric) -- Simulation drag or function handle.
      %     Default: ``0.0``.  Equivalent to no inertia.
      %
      %   - brownian_motion (logical) -- Use brownian motion.
      %     Default: ``true``.
      %
      %   - position (numeric 3x1) -- Example particle position.
      %     Used to estimate drag/mass.  Default: ``[0;0;0]``.
      %
      %   - rotation (numeric 3x3) -- Example particle orientation.
      %     Used to estimate drag/mass.  Default: ``eye(3)``.
      %
      %   - time (numeric) -- Example simulation time.
      %     Used to estimate drag/mass.  Default: ``0.0``.
      %
      %   - temperature (numeric) -- Temperature for Brownian motion
      %     simulations or function handle.  Default: ``300.0`` [Kelvin].
      %
      %   - time_threshold (numeric) -- Time scale threshold used for
      %     method selection.  Default: ``1.0e-6`` [seconds].

      p = inputParser;
      p.addParameter('drag', 1.0);
      p.addParameter('mass', 0.0);
      p.addParameter('brownian_motion', true);
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.addParameter('temperature', 300.0);
      p.addParameter('time_threshold', 1.0e-6);
      p.parse(varargin{:});

      % Get mass and drag
      mass = p.Results.mass;
      drag = p.Results.drag;
      if isa(mass, 'function_handle')
        mass = mass(p.Results.position, p.Results.rotation, p.Results.time);
      end
      if isa(drag, 'function_handle')
        drag = drag(p.Results.position, p.Results.rotation, p.Results.time);
      end

      % Select method using criteria tau = m/gamma < time_threshold
      % This is from https://doi.org/10.1119/1.4772632
      if mass / min(vecnorm(drag, 2, 2)) < p.Results.time_threshold
        sys = ott.dynamics.Stokes(forceMethod, ...
          'drag', p.Results.drag, ...
          'brownian_motion', p.Results.brownian_motion, ...
          'temperature', p.Results.temperature);
      else
        sys = ott.dynamics.Newtonian(forceMethod, ...
          'drag', p.Results.drag, ...
          'mass', p.Results.mass, ...
          'brownian_motion', p.Results.brownian_motion, ...
          'temperature', p.Results.temperature);
      end
    end
  end

  methods
    function sys = DynamicsSystem(forceMethod, varargin)
      % Construct the base class
      %
      % Handles arguments common to both Stokes and Newtonian
      % Perhaps we will remove some of these args/properties in future.

      % Parse arguments
      p = inputParser;
      p.addParameter('drag', 1.0);
      p.addParameter('temperature', 300.0);
      p.addParameter('brownian_motion', true);
      p.parse(varargin{:});

      sys.force_method = forceMethod;

      sys.drag = p.Result.drag;
      sys.temperature = p.results.temperature;
      sys.brownian_motion = brownian_motion;
    end
  end
end

