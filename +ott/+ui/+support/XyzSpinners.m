classdef XyzSpinners < handle
% Creates a labeled set of 3 spinners (e.g., for XYZ coordinates)

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
    XSpinner          matlab.ui.control.Spinner
    YSpinner          matlab.ui.control.Spinner
    ZSpinner          matlab.ui.control.Spinner
  end
  
  methods
    function obj = XyzSpinners(parent, varargin)
      
      p = inputParser();
      p.addParameter('position', [0, 0]);
      p.addParameter('label', 'xyz');
      p.addParameter('visible', 'on');
      p.parse(varargin{:});
      
      sWidth = 60;
      sHeight = 22;
      lWidth = 55;
      
      % Label
      obj.Label = uilabel(parent);
      obj.Label.HorizontalAlignment = 'left';
      obj.Label.Position = [p.Results.position lWidth sHeight];
      obj.Label.Text = p.Results.label;

      % X Spinner
      obj.XSpinner = uispinner(parent);
      obj.XSpinner.Position = [p.Results.position+[lWidth,0], ...
          sWidth sHeight];

      % Y Spinner
      obj.YSpinner = uispinner(parent);
      obj.YSpinner.Position = [p.Results.position+[lWidth+sWidth,0], ...
          sWidth sHeight];

      % Z Spinner
      obj.ZSpinner = uispinner(parent);
      obj.ZSpinner.Position = [p.Results.position+[lWidth+2*sWidth,0], ...
          sWidth sHeight];
        
      % Update visibility
      obj.Visible = p.Results.visible;
    end
  end
  
  methods
    function val = get.Value(obj)
      val = [obj.XSpinner.Value, obj.YSpinner.Value, obj.ZSpinner.Value];
    end
    
    function set.Value(obj, val)
      obj.XSpinner.Value = val(1);
      obj.YSpinner.Value = val(2);
      obj.ZSpinner.Value = val(3);
    end
    
    function val = get.Visible(obj)
      val = obj.XSpinner.Visible;
    end
    
    function set.Visible(obj, val)
      obj.Label.Visible = val;
      obj.XSpinner.Visible = val;
      obj.YSpinner.Visible = val;
      obj.ZSpinner.Visible = val;
    end
  end
end