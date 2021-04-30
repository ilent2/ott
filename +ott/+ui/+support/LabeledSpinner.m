classdef LabeledSpinner < handle
% Creates a labeled spinner

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  events (NotifyAccess = protected)
    ValueChanged
  end

  properties (Dependent)
    Value
    Visible
  end
  
  properties
    Label             matlab.ui.control.Label
    Spinner           matlab.ui.control.Spinner
  end
  
  methods
    function obj = LabeledSpinner(parent, varargin)
      
      p = inputParser();
      p.addParameter('position', [0, 0]);
      p.addParameter('label', 'spinner');
      p.addParameter('visible', 'on');
      p.parse(varargin{:});
      
      sWidth = 130;
      sHeight = 22;
      lWidth = 95;
      
      % Label
      obj.Label = uilabel(parent);
      obj.Label.HorizontalAlignment = 'left';
      obj.Label.Position = [p.Results.position lWidth sHeight];
      obj.Label.Text = p.Results.label;
      obj.Label.Visible = p.Results.visible;

      % Spinner
      obj.Spinner = uispinner(parent);
      obj.Spinner.Position = [p.Results.position+[100,0], ...
          sWidth sHeight];
      obj.Spinner.Visible = p.Results.visible;
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
  end
end