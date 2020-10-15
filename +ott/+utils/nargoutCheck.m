function nargoutCheck(obj, numout)
% Checks the number of arguments is suitable for class property assignment.
%
% Checks if the object is a handle class, if so, has no effect.
% If the object is not a handle, checks numout is at least 1.
%
% Usage
%   ott.utils.nargoutCheck(obj, nargout);

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

if ~isa(obj, 'handle') && numout < 1
  warning('ott:utils:nargoutCheck:no_outputs', ...
      'Expected at least one output argument with non-handle class');
end

