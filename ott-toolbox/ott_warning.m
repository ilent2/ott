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
  internal = false;
end

% Look for ott key words
ott_warning_flag = false;
if length(varargin) == 1
  if strcmpi(varargin{1}, 'once')
    once = true;
    ott_warning_flag = true;
  elseif strcmpi(varargin{1}, 'off')
    once = false;
    internal = false;
  elseif strcmpi(varargin{1}, 'internal')
    internal = true;
    ott_warning_flag = true;
  elseif strcmpi(varargin{1}, 'external')
    internal = false;
    ott_warning_flag = true;
  end
end

if internal
  % Display no warnings
elseif ~ott_warning_flag

  % Call matlab warning function
  varargout{:} = warning(varargin{:});

  % Turn off ott warning if once mode
  if length(varargin) >= 1
    warning_id = varargin{1};
    if once && strcmpi(warning_id(1:3), 'ott')
      warning('off', warning_id);
    end
  end
end
