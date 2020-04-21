classdef Stokes < ott.dynamics.DynamicsSystem
% Identifier and properties for Stokes dynamics (1st order DE, no inertia)
% Inherits from :class:`DynamicsSystem`.
%
% Stores properties for the DE::
%
%   \Gamma \dot{x} = F + F_{BM}
%
% where :math:`\Gamma` is the drag, :math:`\ddot{x}` is the particle
% velocity, :math:`F` is the deterministic forces and :math:`F_{BM}` are
% the non-deterministic forces.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function sys = Stokes(forceMethod, varargin)
      % Create a new identifier for Stokes dynamics
      %
      % Usage
      %   sys = Stokes(forceMethod, ...)
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

      sys = sys@ott.dynamics.DynamicsSystem(forceMethod, varargin{:});
    end
  end
end

