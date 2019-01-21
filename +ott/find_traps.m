function traps = find_traps(position, force, varargin)
% FIND_TRAPS attempt to find and characterise traps from position-force data
%
% traps = find_traps(position, force, ...) attempts to find and characterise
% possible traps for the given position and force data.  Position and
% force should be vectors with position and force along one axis.
%
% The returned traps is an array of structures with information about
% the trap equilibrium position, trap depth and trap stiffness and trap range.
%
% Optional named arguments:
%   keep_unstable    bool   keep unstable equilibriums (default: false)
%   depth_threshold_e num   percentage of max depth for trap acceptance
%       Use [] for no threshold.  (default: 1e-2).
%
% See also ott.find_equilibrium ott.axial_equilibrium and ott.trap_stiffness.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% TODO: We could add options to polyfit both the minmax force and equilibrium

% Parse optional inputs
p = inputParser;
p.addParameter('keep_unstable', false);
p.addParameter('depth_threshold_e', 1e-2);
p.parse(varargin{:});

% This function is not directly concerned with force/torque calculation
ott.warning('ott:findEquilibrium:move', ...
    'This function will move in a future release');

% Check the size of the inputs
assert(isvector(position), 'position must be a vector not a matrix');
assert(isvector(force), 'force must be a vector not a matrix');

% Make sure the vectors are both colum vectors
position = position(:);
force = force(:);

%% Find rough equilibrium positions

% Find first equilibrium position
if force(1) >= 0
  last = find(force < 0, 1);
else
  last = find(force >= 0, 1);
end

% Find remaing equilibriums
eqs = [];
while ~isempty(last)
  eqs(end+1) = last;

  if force(last) >= 0
    last = find(force(last:end) < 0, 1) + last - 1;
  else
    last = find(force(last:end) >= 0, 1) + last - 1;
  end
end

%% Find more precise equilibrium positions and trap stiffness
% This is based on ott.find_equilibrium

peqs = zeros(1, length(eqs));
pstiff = zeros(1, length(eqs));

for ii = 1:length(eqs)

  % Fit polynomial to points aroung equilibrium
  eqrange = max([eqs(ii)-2,1]):min([eqs(ii)+2,length(position)]);
  z = position(eqrange);

  % Scale position and force before polyfit
  zmin = min(z);
  zmax = max(z);
  z = 2 * (z - zmin) / (zmax - zmin) - 1;
  zzero = 2 * (position(eqs(ii)) - zmin) / (zmax - zmin) - 1;

  % Find equilibrium: fit local points to 3rd order polynomail
  % Requires small distance between positions
  pz=polyfit(z, force(eqrange), 3);
  root_z=roots(pz);

  % Ignore non-real roots
  real_z=root_z(imag(root_z)==0);
  if numel(real_z) == 0
    error('No real roots');
  end

  % Keep only one root closest to position(eqs(ii))
  zeqs_idx = abs(real_z-zzero) == min(abs(real_z-zzero));
  real_z = real_z(zeqs_idx);
  real_z = real_z(1);

  % Get the equilibrium
  peqs(ii) = real_z;

  % Calculate stiffness (using derivative of 3rd order polynomial)
  dpz=[3*pz(1),2*pz(2),1*pz(3)];
  pstiff(ii) = polyval(dpz, peqs(ii));

  % Inverse scaling of position and force after fitting
  peqs(ii) = (peqs(ii) + 1)/2*(zmax - zmin) + zmin;
  pstiff(ii) = pstiff(ii)*2/(zmax - zmin);

end

%% Calculate other properties needed for traps

traps = struct('position', {}, 'stiffness', {}, 'depth', {}, ...
    'range', {}, 'minmax_force', {}, 'minmax_position', {});

for ii = 1:length(eqs)

  % Check if stable equilibrium
  if pstiff(ii) >= 0
    continue;
  end
  
  idx = length(traps) + 1;

  % Store stiffness and equilibrium
  traps(idx).position = peqs(ii);
  traps(idx).stiffness = pstiff(ii);

  % Calculate trap range
  traps(idx).range = [-Inf, Inf];
  if ii ~= 1
    traps(idx).range(1) = peqs(ii-1);
  end
  if ii ~= length(eqs)
    traps(idx).range(2) = peqs(ii+1);
  end

  % Calculate trap depth
  frange = [1, length(force)];
  if ii ~= 1
    frange(1) = eqs(ii-1);
  end
  if ii ~= length(eqs)
    frange(2) = eqs(ii+1)-1;
  end
  [fn, fnidx] = min(force(frange(1):frange(2)));
  [fx, fxidx] = max(force(frange(1):frange(2)));
  traps(idx).minmax_force = [fx, fn];
  traps(idx).minmax_position = position(frange(1)+[fxidx, fnidx]-1).';
  traps(idx).depth = min(abs(traps(idx).minmax_force));

end

%% Discard traps that are too shallow

if ~isempty(p.Results.depth_threshold_e)
  max_depth = max(abs(force));
  drop_list = false(size(traps));

  for ii = 1:length(traps)

    if traps(ii).depth < p.Results.depth_threshold_e*max_depth
      drop_list(ii) = true;
    end
  end

  traps(drop_list) = [];
end

