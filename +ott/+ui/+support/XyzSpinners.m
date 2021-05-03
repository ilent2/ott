classdef XyzSpinners < ott.ui.support.LabeledWidget
% Creates a labeled set of 3 spinners (e.g., for XYZ coordinates)

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    Value
    ValueChangedFcn
    Limits
    LowerLimitInclusive
    Step
  end
  
  properties
    XSpinner          matlab.ui.control.Spinner
    YSpinner          matlab.ui.control.Spinner
    ZSpinner          matlab.ui.control.Spinner
  end
  
  methods
    function obj = XyzSpinners(parent, varargin)
      % Create a new XYZ Spinner instance
      
      p = inputParser;
      p.addParameter('label', 'Xyz');
      p.addParameter('visible', 'on');
      p.parse(varargin{:});
      
      obj = obj@ott.ui.support.LabeledWidget(parent, ...
          'label', p.Results.label, 'visible', p.Results.visible);
        
      % Configure grid for additional entries
      obj.Grid.ColumnWidth = {'1x', 60, 60, 60};
      obj.Grid.ColumnSpacing = 0;

      % X Spinner
      obj.XSpinner = uispinner(obj.Grid);
      obj.XSpinner.Layout.Column = 2;
      obj.XSpinner.Layout.Row = 1;

      % Y Spinner
      obj.YSpinner = uispinner(obj.Grid);
      obj.YSpinner.Layout.Column = 3;
      obj.YSpinner.Layout.Row = 1;

      % Z Spinner
      obj.ZSpinner = uispinner(obj.Grid);
      obj.ZSpinner.Layout.Column = 4;
      obj.ZSpinner.Layout.Row = 1;
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
    
    function set.Limits(obj, val)
      obj.XSpinner.Limits = val;
      obj.YSpinner.Limits = val;
      obj.ZSpinner.Limits = val;
    end
    
    function set.Step(obj, val)
      obj.XSpinner.Step = val;
      obj.YSpinner.Step = val;
      obj.ZSpinner.Step = val;
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.XSpinner.ValueChangedFcn = val;
      obj.YSpinner.ValueChangedFcn = val;
      obj.ZSpinner.ValueChangedFcn = val;
    end
    
    function set.LowerLimitInclusive(obj, val)
      obj.XSpinner.LowerLimitInclusive = val;
      obj.YSpinner.LowerLimitInclusive = val;
      obj.ZSpinner.LowerLimitInclusive = val;
    end
  end
end