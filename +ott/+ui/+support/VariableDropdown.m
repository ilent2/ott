classdef VariableDropdown < handle
% A dropdown menu for selecting variables from the current workspace
  
% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Access=public)
    Grid          matlab.ui.container.GridLayout
    Label         matlab.ui.control.Label
    DropDown      matlab.ui.control.DropDown
  end
  
  properties (Dependent)
    Layout
  end
  
  methods
    function obj = VariableDropdown(parent, varargin)
      
      p = inputParser;
      p.parse(varargin{:});
      
      % Create grid
      obj.Grid = uigridlayout(parent, [1, 2]);
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {100, 100};
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
    function val = get.Layout(obj)
      val = obj.Grid.Layout;
    end
    
    function set.Layout(obj, val)
      obj.Grid.Layout = val;
    end
  end
end