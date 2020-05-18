function varargout = RotationPositionDispatcher(func, other, varargin)
% Helper function for processing arrays of objects.
%
% Dispatches calls to a function after translating or rotating
% a particle by a specified position/rotation or an array of
% positions and rotations.
%
% Usage
%   [varargout{1:nargout}] = RotationPositionDispatcher(func, other, ...)
%
% Parameters
%   - func (function_handle) -- Function to call with `other` and
%     any unmatched arguments.
%
%   - other -- An object with rotation and position properties.
%
% Named arguments
%   - position (3xN numeric) -- Position to apply to other before
%     dispatching to func.  Default: ``[]``.
%
%   - rotation (3x3N numeric) -- Rotation to apply to other before
%     dispatching to func.  Default: ``[]``.
%
%   - prxcat (numeric|empty) -- Dimension to concatenate prxfun
%     results along.  Default: ``[]``.
%
% Length of position and rotation must match or be equal to 1.

% TODO: Should this be its own function or merge with prxfun?
% TODO: This might move again soon

assert(isa(other, 'ott.utils.RotationPositionEntity'), ...
    'other must be a rotationpositionentity');

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('position', []);
p.addParameter('rotation', []);
p.addParameter('prxcat', []);
p.parse(varargin{:});
unmatched = ott.utils.unmatchedArgs(p);

% Check if we need to operate on multiple beams
if size(p.Results.position, 2) > 1 || size(p.Results.rotation, 2) > 3
  [varargout{1:nargout}] = ott.utils.prxfun(@(varargin) ...
      ott.utils.RotationPositionProp.rotationPositionHelper(...
      func, other, varargin{:}), ...
      'position', p.Results.position, ...
      'rotation', p.Results.rotation, unmatched{:});

  % Concatenate outputs along requested dimension
  if ~isempty(p.Results.prxcat)
    for ii = 1:nargout
      varargout{ii} = cat(p.Results.prxcat, varargout{ii}{:});
    end
  end

  return;
end

% Add position/rotation to beam
if ~isempty(p.Results.position)
  other.position = other.position - p.Results.position;
end
if ~isempty(p.Results.rotation)
  other.rotation = p.Results.rotation.' * other.rotation;
end

[varargout{1:nargout}] = func(other, unmatched{:});

