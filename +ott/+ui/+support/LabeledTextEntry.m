classdef LabeledTextEntry < ott.ui.support.LabeledWidget
% A labeled text entry widget

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    Value
    ValueChangedFcn
  end
  
  properties
    Entry           matlab.ui.control.EditField
  end
  
  methods
    function obj = LabeledTextEntry(parent, varargin)
      % Create a labeled spinner widget
      
      p = inputParser;
      p.addParameter('label', 'Entry');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      obj = obj@ott.ui.support.LabeledWidget(parent, ...
        'label', p.Results.label, unmatched{:});

      % Spinner
      obj.Entry = uieditfield(obj.Grid, 'text');
      obj.Entry.Layout.Column = 2;
      obj.Entry.Layout.Row = 1;
      
    end
  end
  
  methods
    function val = get.Value(obj)
      val = obj.Entry.Value;
    end
    
    function set.Value(obj, val)
      obj.Entry.Value = val;
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.Entry.ValueChangedFcn = val;
    end
  end
end