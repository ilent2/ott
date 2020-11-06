function m = arrayMaxSize()
% Get the maximum array size (in doubles) that fits into physical memory.
%
% Internally, uses ``memory`` to get the available physical memory.
% On some systems this fails, in this case, we simply raise
% a warning 'ott:utils:arrayMaxSize:cant_get_memory' and return ``Inf``.
% This is only a guide to available memory and doesn't include swap.
%
% Usage
%   m = ott.utils.arrayMaxSize()

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

try
  [uV, sV] = memory;
catch
  m = Inf;
  warning('ott:utils:arrayMaxSize:cant_get_memory', ...
      'Unable to get array max size using memory command');
  return;
end

m = floor(sV.PhysicalMemory.Available/8);

