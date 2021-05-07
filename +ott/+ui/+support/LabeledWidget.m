classdef (Abstract) LabeledWidget < ott.ui.support.GridWidget
% Base class for labeled widgets
%
% Abstract properties
%   - Value
%   - ValueChangedFcn
%
% Dependent properties
%   - Visible       -- Shows/hides the containing grid
%   - Enable        -- Sets the enable property of all child widgets

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.
  
  properties (Dependent)
    Visible
    Enable
  end
  
  properties (Abstract)
    Value
    ValueChangedFcn
  end
  
  properties
    Label             matlab.ui.control.Label
  end
  
  methods
    function obj = LabeledWidget(parent, varargin)
      % Create a new labeled widget template
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      p = inputParser();
      p.addParameter('label', 'Label');
      p.addParameter('visible', 'on');
      p.addParameter('wwidth', 100);
      p.parse(varargin{:});
      
      % Configure grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {'1x', p.Results.wwidth};
      obj.Grid.ColumnSpacing = 1;
      obj.Grid.RowSpacing = 1;
      
      % Label
      obj.Label = uilabel(obj.Grid);
      obj.Label.HorizontalAlignment = 'left';
      obj.Label.Text = p.Results.label;
      obj.Label.Layout.Column = 1;
      obj.Label.Layout.Row = 1;
      
      % Set visibility
      obj.Visible = p.Results.visible;
      
    end
  end
  
  methods
    function val = get.Visible(obj)
      val = obj.Grid.Visible;
    end
    
    function set.Visible(obj, val)
      obj.Grid.Visible = val;
    end
    
    function set.Enable(obj, val)
      widgets = findobj(obj.Grid.Children, '-property', 'Enable');
      set(widgets, 'Enable', val);
    end
  end
end