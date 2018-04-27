function change_warnings(enable)
%CHANGE_WARNINGS enables or disables move/depreciation warnings
%
% CHANGE_WARNINGS 'on' or CHANGE_WARNINGS('on') enables warnings.
%
% CHANGE_WARNINGS 'off' or CHANGE_WARNINGS('off') disables warnings.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Check number of inputs
if nargin == 0
  error('ott:change_warnings:no_input', ...
      'change_warnings expects input: on or off to enable/disable warnings');
end

% Check input is on or off
if ~strcmpi(enable, 'on') && ~strcmpi(enable, 'off')
  error('ott:change_warnings:input_value', ...
      'change_warnings expects input to be either on or off');
end

% List of warnings we control
warnings = { 'ott:findEquilibrium:move', ...
    'ott:axialEquilibrium:move' };

% Turn off warnings
for ii = 1:length(warnings)
  warning(enable, warnings{ii});
end
