classdef Newtonian < ott.dynamics.DynamicsSystem
% Identifier and properties for Newtonian dynamics (2nd order DE, inertia)
% Inherits from :class:`DynamicsSystem`.
%
% Stores properties for the DE::
%
%   m \ddot{x} = -\Gamma \dot{x} + F + F_{BM}
%
% where :math:`m` is the particle mass, :math:`\Gamma` is the drag,
% :math:`\ddot{x}` is the particle velocity, :math:`F` is the
% deterministic forces and :math:`F_{BM}` are the non-deterministic forces.
%
% Properties
%   - force_method          -- Force calculation method
%   - drag                  -- Drag tensor
%   - temperature           -- Temperature for brownian motion
%   - brownian_motion       -- If brownian motion is used
%   - mass                  -- Particle mass

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    mass                    % Particle mass
  end

  methods
    function sys = Newtonian(forceMethod, varargin)
      % Create a new identifier for Newtonian dynamics
      %
      % Usage
      %   sys = Newtonian(forceMethod, ...)
      %
      % Optional named arguments
      %   - drag -- Drag or function handle
      %     Default: ``1.0``.
      %
      %   - temperature -- Temperature or function handle
      %     Default: ``300.0`` [Kelvin].
      %
      %   - brownian_motion (logical) - If brownian motion should be used.
      %     Default: ``true``.
      %
      %   - mass -- Mass or function handle.
      %     Default: ``1.0``.

      p = inputParser;
      p.keepUnmatched = true;
      p.addParameter('mass', 1.0);
      p.parse(varargin{:});

      % Construct base
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];
      sys = sys@ott.dynamics.DynamicsSystem(forceMethod, unmatched{:});

      % Add mass
      sys.mass = p.Results.mass;
    end
  end
end

