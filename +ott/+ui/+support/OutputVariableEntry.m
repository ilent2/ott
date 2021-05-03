classdef OutputVariableEntry < ott.ui.support.GridWidget
% A Label and Text Entry for specifying the output variable name.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    Label         matlab.ui.control.Label
    EditField     matlab.ui.control.EditField
  end
  
  properties (Dependent)
    Value
    ValueChangedFcn
  end
  
  methods
    function obj = OutputVariableEntry(parent, varargin)
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      p = inputParser;
      p.parse(varargin{:});
      
      % Create grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {'1x', 100};
      obj.Grid.ColumnSpacing = 1;
      obj.Grid.RowSpacing = 1;
      
      % Create label
      obj.Label = uilabel(obj.Grid);
      obj.Label.HorizontalAlignment = 'left';
      obj.Label.Layout.Column = 1;
      obj.Label.Layout.Row = 1;
      obj.Label.Text = 'Variable Name';

      % Create ScatteredBeamEditField
      obj.EditField = uieditfield(obj.Grid, 'text');
      obj.EditField.Layout.Column = 2;
      obj.EditField.Layout.Row = 1;
      obj.EditField.Value = '';
      
    end
    
    function WriteVariable(obj, data)
      % Write the specified data to the given variable name
      
      % Check we have work to do
      if isempty(obj.Value)
        error('No variable name specified');
      end
      
      % Check valid variable name
      if ~isvarname(obj.Value)
        error('Value must be a valid variable name');
      end
      
      % Store result
      assignin('base', obj.Value, data);
      
    end
  end
  
  methods
    function val = get.Value(obj)
      val = obj.EditField.Value;
    end
    
    function set.Value(obj, val)
      obj.EditField.Value = val;
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.EditField.ValueChangedFcn = val;
    end
  end
end