classdef RangeSpinners < ott.ui.support.LabeledWidget
% Creates spinners for inputting a range.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Value checking

  properties (Dependent)
    Value
    ValueChangedFcn
    Limits
    Step
  end
  
  properties
    LowerSpinner    matlab.ui.control.Spinner
    UpperSpinner    matlab.ui.control.Spinner
  end
  
  methods
    function obj = RangeSpinners(parent, varargin)
      % Create a new range spinner instance
      
      p = inputParser;
      p.addParameter('label', 'Range');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      obj = obj@ott.ui.support.LabeledWidget(parent, ...
        'label', p.Results.label, unmatched{:});
        
      % Configure grid for additional entries
      wwidth = obj.Grid.ColumnWidth{2};
      obj.Grid.ColumnWidth = {'1x', wwidth/2, wwidth/2};
      obj.Grid.ColumnSpacing = 0;

      % X Spinner
      obj.LowerSpinner = uispinner(obj.Grid);
      obj.LowerSpinner.Layout.Column = 2;
      obj.LowerSpinner.Layout.Row = 1;

      % Y Spinner
      obj.UpperSpinner = uispinner(obj.Grid);
      obj.UpperSpinner.Layout.Column = 3;
      obj.UpperSpinner.Layout.Row = 1;
      
    end
  end
  
  methods % Getters/setters
    function val = get.Value(obj)
      val = [obj.LowerSpinner.Value, obj.UpperSpinner.Value];
    end
    
    function set.Value(obj, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'Value must be 2 element numeric');
      obj.LowerSpinner.Value = val(1);
      obj.UpperSpinner.Value = val(2);
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.LowerSpinner.ValueChangedFcn = val;
      obj.UpperSpinner.ValueChangedFcn = val;
    end
    
    function set.Limits(obj, val)
      obj.LowerSpinner.Limits = val;
      obj.UpperSpinner.Limits = val;
    end
    
    function set.Step(obj, val)
      obj.LowerSpinner.Step = val;
      obj.UpperSpinner.Step = val;
    end
  end
end