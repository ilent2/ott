function varargout = prxfun(func, sz, varargin)
% Position-rotation expansion function
%
% Applies a specified function to each combination of the input
% positions and rotations.
%
% Usage
%   ret = prxfun(func, sz, ...)
%
%   [ret1, ret2, ...] = prxfun(func, sz, ...)
%   Returns multiple arguments.  Requires func to return multiple results.
%
% Parameters
%   - func (function_handle) -- The function to call for each combination
%     of the position and rotation parameters.  The function should
%     take named arguments for position and rotation.
%
%   - sz (numeric) -- Size of the output returned by `func`.
%
% Optional named arguments
%   - position (3xN numeric) -- position arguments.  Default: ``[]``.
%   - rotation (3x3N numeric) -- rotation arguments. Default: ``[]``.
%   - progress (function_handle) -- a function to call after each step.
%     Default: ``[]``.  The function handle should take two arguments,
%     the current step index and the total amount of work.
%
%   - zeros (function_handle) -- A function like ``zeros`` to call to
%     allocate memory for the output.  Default: ``zeros``.
%
% Unmatched arguments are passed to func.
% Number of positions and rotations must be scalar or equal.

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('position', []);
p.addParameter('rotation', []);
p.addParameter('progress', []);
p.addParameter('zeros', @zeros);
p.parse(varargin{:});
unmatched = ott.utils.unmatchedArgs(p);

position = p.Results.position;
rotation = p.Results.rotation;

% Check input types
assert(isnumeric(position), 'Position must be numeric');
assert(isnumeric(rotation), 'Rotation must be numeric');

% Check size of inputs
assert(size(position, 1) == 0 || size(position, 1) == 3, ...
    'Position must be empty or be 3xN matrix');
assert((size(rotation, 1) == 0 || size(rotation, 1) == 3) ...
    && mod(size(rotation, 2), 3) == 0, ...
    'Position must be empty or be 3x3N matrix');

% Get amount of work
Nposition = size(position, 2);
Nrotation = size(rotation, 2)/3;
Nwork = max(Nposition, Nrotation);

% Check size of work matches
assert(Nposition == Nrotation || Nposition == 0 || Nrotation == 0 ...
    || Nposition == 1 || Nrotation == 1, ...
    'Length of position/rotation must be empty, scalar or equal');

% If there is no work, just call the function
if Nwork == 0
  [varargout{1:nargout}] = func('position', [], 'rotation', []);
  return
end

% Allocate memory for output(s)
if iscell(p.Results.zeros)
  for ii = 1:length(p.Results.zeros)
    varargout{ii} = p.Results.zeros{ii}([sz, Nwork]);
  end
else
  varargout = repmat({p.Results.zeros([sz, Nwork])}, 1, nargout);
end

epw = prod(sz);

for ii = 1:Nwork

  % Get the position
  if Nposition == 1 || Nposition == 0
    our_position = position;
  else
    our_position = position(:, ii);
  end

  % Get the rotation
  if Nrotation == 1 || Nrotation == 0
    our_rotation = rotation;
  else
    our_rotation = rotation(:, (1:3) + (ii-1)*3);
  end

  % Calculate the work
  func_outputs = {};
  [func_outputs{1:nargout}] = func('position', our_position, ...
      'rotation', our_rotation, unmatched{:});
  for jj = 1:length(func_outputs)
    varargout{jj}((1:epw) + (ii-1)*epw) = func_outputs{jj};
  end

  if ~isempty(p.Results.progress)
    p.Results.progress(ii, Nwork);
  end
end


