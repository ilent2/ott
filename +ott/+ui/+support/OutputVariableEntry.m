classdef OutputVariableEntry < handle
% A Label and Text Entry for specifying the output variable name.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    Grid          matlab.ui.container.GridLayout
    Label         matlab.ui.control.Label
    EditField     matlab.ui.control.EditField
  end
  
  properties (Dependent)
    Layout
  end
  
  methods
    function obj = OutputVariableEntry(parent, varargin)
      
      p = inputParser;
      p.parse(varargin{:});
      
      % Create grid
      obj.Grid = uigridlayout(parent, [1, 2]);
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {100, 100};
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
  end
  
  methods
    function val = get.Layout(obj)
      val = obj.Grid.Layout;
    end
    
    function set.Layout(obj, val)
      obj.Grid.Layout = val;
    end
  end
end