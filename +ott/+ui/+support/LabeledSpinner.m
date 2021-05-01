classdef LabeledSpinner < ott.ui.support.GridWidget
% Creates a labeled spinner

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  events (NotifyAccess = protected)
    ValueChanged
  end

  properties (Dependent)
    Value
    ValueChangedFcn
    Visible
    Step
    Limits
  end
  
  properties
    Label             matlab.ui.control.Label
    Spinner           matlab.ui.control.Spinner
  end
  
  methods
    function obj = LabeledSpinner(parent, varargin)
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      p = inputParser();
      p.addParameter('label', 'spinner');
      p.addParameter('visible', 'on');
      p.parse(varargin{:});
      
      % Create grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {'1x', 100};
      obj.Grid.ColumnSpacing = 1;
      obj.Grid.RowSpacing = 1;
      
      % Label
      obj.Label = uilabel(obj.Grid);
      obj.Label.HorizontalAlignment = 'left';
      obj.Label.Text = p.Results.label;
      obj.Label.Layout.Column = 1;
      obj.Label.Layout.Row = 1;

      % Spinner
      obj.Spinner = uispinner(obj.Grid);
      obj.Spinner.Layout.Column = 2;
      obj.Spinner.Layout.Row = 1;
      
      % Set visibility
      obj.Visible = p.Results.visible;
    end
  end
  
  methods
    function val = get.Value(obj)
      val = obj.Spinner.Value;
    end
    
    function set.Value(obj, val)
      obj.Spinner.Value = val;
    end
    
    function val = get.Visible(obj)
      val = obj.Spinner.Visible;
    end
    
    function set.Visible(obj, val)
      obj.Spinner.Visible = val;
      obj.Label.Visible = val;
    end
    
    function set.Limits(obj, val)
      obj.Spinner.Limits = val;
    end
    
    function set.Step(obj, val)
      obj.Spinner.Step = val;
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.Spinner.ValueChangedFcn = val;
    end
  end
end