classdef VariableDropdown < ott.ui.support.LabeledDropDown
% A dropdown menu for selecting variables from the current workspace
  
% Copyright 2020 Isaac Lenon.
% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: We could look at somthing like this for input validation:
%   https://undocumentedmatlab.com/articles/editbox-data-input-validation

  properties (Access=public)
    Filter (1, :) char
  end
  
  properties (Dependent, SetAccess=private)
    Variable
  end
  
  methods
    function obj = VariableDropdown(parent, varargin)
      
      p = inputParser;
      p.addParameter('label', 'Variable');
      p.addParameter('filter', 'double');
      p.addParameter('visible', 'on');
      p.parse(varargin{:});
      
      obj = obj@ott.ui.support.LabeledDropDown(parent, ...
        'visible', p.Results.visible, 'label', p.Results.label);
      
      obj.Filter = p.Results.filter;

      % Configure DropDown
      obj.DropDown.Editable = 'on';
      obj.DropDown.BackgroundColor = [1 1 1];
      obj.DropDown.Value = {};
      
    end
    
    function refresh(obj)
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

      if nargin < 3
        basevars = {};
      end

      varnames = evalin('base', 'who');
      vars = basevars;
      for ii = 1:length(varnames)
          if isa(evalin('base', varnames{ii}), obj.Filter)
              vars{end+1} = varnames{ii}; %#ok<AGROW>
          end
      end

      % If present value not found, unset the selection
      if ~isempty(obj.DropDown.Value) && ~any(strcmpi(obj.DropDown.Value, vars))
        obj.DropDown.Value = '';
      end

      % Populate DropDown with device list
      obj.DropDown.Items = vars;
      obj.DropDown.ItemsData = vars;
    end
  end
  
  methods
    function val = get.Variable(obj)
      try
        val = evalin('base', obj.Value);
      catch
        val = [];
      end
    end
  end
end