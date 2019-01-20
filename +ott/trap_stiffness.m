function [kf, kt, K] = trap_stiffness(beam, T, varargin)
% TRAP_STIFFNESS calculate the force and torque trap stiffness
%
% [kf, kt, K] = TRAP_STIFFNESS(beam, T) calculate the trap stiffness
% for the force, kf, and torque, kt for a particle, T, in a beam.
% K is the full 6x2N stiffness matrix, each column corresponds to a
% different direction.  K can be multiplied by the drag tensor to
% calculate the fluid+optical trap stiffness.
%
% TRAP_STIFFNESS(..., 'position', x0) and TRAP_STIFFNESS(..., 'rotation', R)
% specify the particle position and rotation.
%
% TRAP_STIFFNESS(..., 'direction', d) specifies the directions to
% calculate the trap stiffness in.  Specify 'beam' for the XYZ direction
% for the current beam orientation.  For specific directions, specify
% a 3xN matrix of vectors to calculate the force along and torque around.
% Default is 'beam'.
%
% TRAP_STIFFNESS(..., 'method', m) specifies the method to use for
% calculating the trap stiffness.
% Supported methods [calculations direction/other]:
%     'cntr'    use central differences [4/0]
%     'fwd'     use forward differences [2/1]
%
% TRAP_STIFFNESS(..., 'step', [ dx dt ]) specifies the position and
% torque calculation step sizes (in beam units and radians).
% Default: 1e-3/beam.k_medium and 1e-3.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

p = inputParser;
p.addParameter('direction', 'beam');
p.addParameter('method', 'cntr');
p.addParameter('rotation', eye(3));
p.addParameter('position', [0;0;0]);
p.addParameter('step', [ 1e-3/beam.k_medium 1e-3 ]);
p.parse(varargin{:});

% Rotate the beam
if ~isempty(p.Results.rotation)
  beam = beam.rotate(p.Results.rotation);
end

% Get the directions to translate/rotated about
D = p.Results.direction;
if strcmpi(D, 'beam')
  D = p.Results.rotation' * eye(3);
elseif strmpi(D, 'particle')
  D = eye(3);
end

% Assign names to step parameters
dx = p.Results.step(1);
dt = p.Results.step(2);

% Allocate space for the output
K = zeros(6, 2*size(D, 2));

if strcmpi(p.Results.method, 'cntr')

  % Create matrix of positions
  x = [ D * -dx, D * dx ] + p.Results.position;

  % Create matrix of all rotations
  R = zeros(3, 3*size(D, 2)*2);
  for ii = 1:size(D, 2)
    R(:, (1:3) + 3*(ii-1)) = ott.utils.rotation_matrix(D(:, ii), dt);
    R(:, (1:3) + 3*(ii-1 + size(D, 2))) = R(:, (1:3) + 3*(ii-1))';
  end

  idx1 = 1:size(D, 2);
  idx2 = size(D, 2)+idx1;

  % Calculate the force stiffness
  [ff, tt] = ott.forcetorque(beam, T, 'position', x);
  K(1:3, idx1) = (ff(:, idx2) - ff(:, idx1))/2/dx;
  K(4:6, idx1) = (tt(:, idx2) - tt(:, idx1))/2/dx;

  % Calculate the rotation about this axis
  [ff, tt] = ott.forcetorque(beam, T, 'rotation', R, ...
      'position', p.Results.position);
  K(1:3, idx2) = (ff(:, idx2) - ff(:, idx1))/2/dt;
  K(4:6, idx2) = (tt(:, idx2) - tt(:, idx1))/2/dt;

elseif strcmpi(p.Results.method, 'fwd')

  % Calculate the common force/torque
  [ff0, tt0] = ott.forcetorque(beam, T, 'position', p.Results.position);

  % Create matrix of all rotations
  R = zeros(3, 3*size(D, 2));
  for ii = 1:size(D, 2)
    R(:, (1:3) + 3*(ii-1)) = ott.utils.rotation_matrix(D(:, ii), dt);
  end

  idx1 = 1:size(D, 2);
  idx2 = size(D, 2)+idx1;

  % Calculate the force stiffness
  [ff, tt] = ott.forcetorque(beam, T, ...
      'position', D*dx + p.Results.position);
  K(1:3, idx1) = (ff - ff0)/dx;
  K(4:6, idx1) = (tt - tt0)/dx;

  % Calculate the rotation about this axis
  [ff, tt] = ott.forcetorque(beam, T, 'rotation', R, ...
      'position', p.Results.position);
  K(1:3, idx2) = (ff - ff0)/dt;
  K(4:6, idx2) = (tt - tt0)/dt;

else
  error('Unsupported method');
end

% Calculate kf and kt
kf = sum(K(1:3, 1:size(D, 2)) .* D, 1).';
kt = sum(K(4:6, size(D, 2)+1:end) .* D, 1).';
