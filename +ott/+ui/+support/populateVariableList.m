function populateVariableList(list, type_name, basevars)
% Populates the variable list with variables of a certain type
%
% Usage
%   populateVariableList(list, type_name) populates the drop down list
%   with variables from the base workspace matching the given type.
%
% Parameters
%   - list (uidropdown) -- List handle to add items to
%   - type_name (type name) -- Name of type to filter variables.
%
% This function is used to populate the contents of a ``uidropdown``
% widget. The function takes a handle to the ``uidropdown`` widget, and
% a Matlab class name and searches the base workspace for variables
% with the specified type. For example usage,
% see :class:`ui.tools.ForceProfile`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  if nargin < 3
    basevars = {};
  end

  varnames = evalin('base', 'who');
  vars = basevars;
  for ii = 1:length(varnames)
      if isa(evalin('base', varnames{ii}), type_name)
          vars{end+1} = varnames{ii}; %#ok<AGROW>
      end
  end
  
  % If present value not found, unset the selection
  if ~isempty(list.Value) && ~any(strcmpi(list.Value, vars))
    list.Value = '';
  end
  
  % Populate DropDown with device list
  list.Items = vars;
  list.ItemsData = vars;
  
end
