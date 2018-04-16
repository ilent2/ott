function varargout = ott_warning(varargin)
%OTT_WARNING overload of MATLAB warning function for OTT
%
% Adds additional options for warnings such as displaying warnings
% once per session and suppressing internal warnings.
%
% OTT_WARNING('internal') and OTT_WARNING('external') are for use by
% ott functions that call other function which likely raise warnings.
%
% OTT_WARNING('once') makes warnings only display once before the warning
% is turned off.  Turning off all warnings resets the once flag.
%
% All other commands are passed onto the builtin warning function.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Setup persistent variables for tracking state
persistent once internal;
if isempty(once)
  once = false;
  internal = 0;
end

% Look for ott key words
ott_warning_flag = false;
if length(varargin) == 1
  if strcmpi(varargin{1}, 'once')
    once = true;
    ott_warning_flag = true;
  elseif strcmpi(varargin{1}, 'off')
    once = false;
    internal = 0;
  elseif strcmpi(varargin{1}, 'internal')
    internal = internal + 1;
    ott_warning_flag = true;
  elseif strcmpi(varargin{1}, 'external')
    internal = max(internal - 1, 0);
    ott_warning_flag = true;
  end
end

if internal ~= 0
  % Display no warnings
elseif ~ott_warning_flag

  % Call matlab warning function
  if nargout ~= 0
    varargout{:} = warning(varargin{:});
  else
    % This case is for octave
    warning(varargin{:});
  end

  % Turn off ott warning if once mode
  if length(varargin) >= 1
    warning_id = varargin{1};
    if once && length(warning_id) > 3 && strcmpi(warning_id(1:3), 'ott')
      warning('off', warning_id);
    end
  end
end
