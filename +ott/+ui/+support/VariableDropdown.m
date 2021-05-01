classdef VariableDropdown < ott.ui.support.GridWidget
% A dropdown menu for selecting variables from the current workspace
  
% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Access=public)
    Label         matlab.ui.control.Label
    DropDown      matlab.ui.control.DropDown
  end
  
  properties (Dependent)
    Value
    ValueChangedFcn
  end
  
  methods
    function obj = VariableDropdown(parent, varargin)
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      p = inputParser;
      p.parse(varargin{:});
      
      % Configure grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {'1x', 100};
      obj.Grid.ColumnSpacing = 1;
      obj.Grid.RowSpacing = 1;
      
      % Create BeamDropDownLabel
      obj.Label = uilabel(obj.Grid);
      obj.Label.Layout.Column = 1;
      obj.Label.Layout.Row = 1;
      obj.Label.Text = 'Variable';

      % Create BeamDropDown
      obj.DropDown = uidropdown(obj.Grid);
      obj.DropDown.Items = {};
      obj.DropDown.Editable = 'on';
      obj.DropDown.BackgroundColor = [1 1 1];
      obj.DropDown.Layout.Column = 2;
      obj.DropDown.Layout.Row = 1;
      obj.DropDown.Value = {};
      
    end
  end
  
  methods
    function set.Value(obj, val)
      obj.DropDown.Value = val;
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.DropDown.ValueChangedFcn = val;
    end
  end
end