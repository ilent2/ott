function R = euler2rot(varargin)
% Convert an euler rotation (rotx, roty, rotz) to a rotation matrix.
%
% This applies rotations in the order: Rz*Ry*Rx*vec, i.e., the same
% form used in :class:`RotateHelper`.
%
% Usage
%   rot = euler2rot(theta_x, theta_y, theta_z, ...)
%
%   rot = euler2rot([theta_x; theta_y; theta_z], ...)
%
% Optional named arguments:
%   usecell    bool     True to output as cell array instead of 3xN matrix.
%       Default: false.  The cell array has the same shape as angle_deg.

% Copyright 2020 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

p = inputParser;
p.addOptional('x', [], @isnumeric);
p.addOptional('y', [], @isnumeric);
p.addOptional('z', [], @isnumeric);
p.addParameter('usecell', false);
p.parse(varargin{:});

if ~any(strcmpi(p.UsingDefaults, 'x')) ....
    && any(strcmpi(p.UsingDefaults, 'y')) ...
    && any(strcmpi(p.UsingDefaults, 'z'))
  
  euler = p.Results.x;
  assert(ismatrix(euler) && isnumeric(euler) && size(euler, 1) == 3, ...
    'Input must be a 3xN numeric matrix.');
  
  Rx = ott.utils.rotx(euler(1, :), 'usecell', true);
  Ry = ott.utils.roty(euler(2, :), 'usecell', true);
  Rz = ott.utils.rotz(euler(3, :), 'usecell', true);
  
elseif ~any(strcmpi(p.UsingDefaults, 'x')) ....
    && ~any(strcmpi(p.UsingDefaults, 'y')) ...
    && ~any(strcmpi(p.UsingDefaults, 'z'))
  theta_x = p.Results.x;
  theta_y = p.Results.y;
  theta_z = p.Results.z;
  
  assert(isnumeric(theta_y), ...
    'Input x  must be numeric');
  assert(isnumeric(theta_y) && numel(theta_x) == numel(theta_y), ...
    'Input y  must be numeric and number of elements must match input x');
  assert(isnumeric(theta_z) && numel(theta_x) == numel(theta_z), ...
    'Input z  must be numeric and number of elements must match input x');
  
  Rx = ott.utils.rotx(theta_x, 'usecell', true);
  Ry = ott.utils.roty(theta_y, 'usecell', true);
  Rz = ott.utils.rotz(theta_z, 'usecell', true);
else
  error('Expected 1 or 3 input arguments');
end

R = cellfun(@(x,y,z) x*y*z, Rx, Ry, Rz, 'UniformOutput', false);

if ~p.Results.usecell
  R = cell2mat(R(:).');
end
