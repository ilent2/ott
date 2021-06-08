classdef LabeledSpinner < ott.ui.support.LabeledWidget
% Creates a labeled spinner

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    Value
    ValueChangedFcn
    Step
    Limits
    LowerLimitInclusive
    RoundFractionalValues
  end
  
  properties
    Spinner           matlab.ui.control.Spinner
  end
  
  methods
    function obj = LabeledSpinner(parent, varargin)
      % Create a labeled spinner widget
      
      p = inputParser;
      p.addParameter('label', 'Spinner');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      obj = obj@ott.ui.support.LabeledWidget(parent, ...
        'label', p.Results.label, unmatched{:});

      % Spinner
      obj.Spinner = uispinner(obj.Grid);
      obj.Spinner.Layout.Column = 2;
      obj.Spinner.Layout.Row = 1;
      
    end
  end
  
  methods
    function val = get.Value(obj)
      val = obj.Spinner.Value;
    end
    
    function set.Value(obj, val)
      obj.Spinner.Value = val;
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
    
    function set.LowerLimitInclusive(obj, val)
      obj.Spinner.LowerLimitInclusive = val;
    end
    
    function set.RoundFractionalValues(obj, val)
      obj.Spinner.RoundFractionalValues = val;
    end
  end
end