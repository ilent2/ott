function varargout = position_force_curves(beam, tmatrix, varargin)
% Calculate position-force curves for a particle in a beam.
%
% The curve is calculated by displacing the particle by a fixed amount
% and calculating the force.
% Not to be confused with :func:`force_position_curves`.
%
% Usage:
%   position_force_curves(beam, tmatrix, ...) plots the position-force
%   curves in the current figure handle.
%
%   [x, fx, y, fy, z, fz] = position_force_curves(beam, tmatrix, ...)
%   Calculates the position-force curves and returns the positions and
%   corresponding forces along each dimension.
%   Outputs 6 3xN matrices for the position and force at each position.
%   If `calculate_torque` is true, fx/fy/fz are 6xN matrices, the last
%   three rows are the calculated torque.
%
% Optional named arguments:
%     - calculate_torque  (logical) -- Calculate (and show) the torque.
%       Default: `false`.
%     - show_plots (logical)        -- Draw the plots to the current figure.
%       Default: `nargout == 0`.
%
% Example::
%   beam = ott.bsc.PmGauss();
%   tmatrix = ott.tmatrix.Mie();
%   ott.tools.position_force_curves(beam, tmatrix);
%
% See also ott.tools.find_traps, ott.bsc.PmGauss, ott.tmatrix.Mie.

p = inputParser();
p.addParameter('origin', []);
p.addParameter('orientation', []);
p.addParameter('calculate_torque', false);
p.addParameter('show_plots', nargout == 0);
p.parse(varargin{:});

% TODO: calculate_torque and show_plots
% TODO: Use orientation and origin
% TODO: Axis orientation vs particle rotation (relative to beam?)
% TODO: Option for only restoring force in plot

% Calculate position locations
x = linspace(-1, 1, 100)*beam.wavelength;
y = linspace(-1, 1, 100)*beam.wavelength;
z = linspace(-1, 1, 100)*beam.wavelength;

% Convert to 3-vectors
x = [x; 0*x; 0*x];
y = [0*y; y; 0*y];
z = [0*z; 0*z; z];

% Calculate forces
fx = ott.forcetorque(beam, tmatrix, 'position', x);
fy = ott.forcetorque(beam, tmatrix, 'position', y);
fz = ott.forcetorque(beam, tmatrix, 'position', z);

% TODO: What about torque?

% Options for axies to include

% TODO

% If beam has a maximum displacement: i.e. plane wave types, only
% calculate force up to this maximum.

% Otherwise, calculate force until force vanishes, or until some
% upper limit is reached.

% Calculate should be efficient, only calculate points where the
% force-curves vary rapidly.  We should be able to verify this (in
% a unit test) by looking at the interpolated force curve.
