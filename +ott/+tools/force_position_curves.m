function varargout = force_position_curves(beam, tmatrix, varargin)
% Calculate particle position after applying a constant force.
%
% This function applies series of constant forces to the particle
% and after allowing the particle to equilibrate, records the position.
% Useful for calculating escape trajectories.
% Not to be confused with :func:`position_force_curves`.
%
% Usage:
%   force_position_curves(beam, tmatrix, ...) calculates the steady-state
%   escape trajectory along the x-axis and plots it in the current figure.
%
%   [f, x] = force_position_curves(beam, tmatrix, ...) calculates the
%   steady-state escape trajectory and returns the force and position values.
%
% Optional named parameters:
%   - direction (3-vector)   -- direction to apply force.
%     Default: `[1;0;0]`.
%   - max_eq_time (scalar)   -- maximum time to wait for the particle
%     to reach equilibrium.
%     Default: `Inf`.
%
% Example::
%   %   beam = ott.bsc.PmGauss();
%   tmatrix = ott.tmatrix.Mie();
%   ott.tools.force_position_curves(beam, tmatrix);
%
% See also ott.tools.find_traps, ott.bsc.PmGauss, ott.tmatrix.Mie.

% TODO

