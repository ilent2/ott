classdef LabeledDropDown < ott.ui.support.LabeledWidget
% A labeled dropdown widget

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    Value
    ValueChangedFcn
    Items
    ItemsData
  end
  
  properties
    DropDown           matlab.ui.control.DropDown
  end
  
  methods
    function obj = LabeledDropDown(parent, varargin)
      % Create a labeled spinner widget
      
      p = inputParser;
      p.addParameter('label', 'DropDown');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      obj = obj@ott.ui.support.LabeledWidget(parent, ...
        'label', p.Results.label, unmatched{:});

      % Spinner
      obj.DropDown = uidropdown(obj.Grid);
      obj.DropDown.Layout.Column = 2;
      obj.DropDown.Layout.Row = 1;
      
    end
  end
  
  methods
    function val = get.Value(obj)
      val = obj.DropDown.Value;
    end
    
    function set.Value(obj, val)
      obj.DropDown.Value = val;
    end
    
    function set.Items(obj, val)
      obj.DropDown.Items = val;
    end
    
    function val = get.ItemsData(obj)
      val = obj.DropDown.ItemsData;
    end
    
    function set.ItemsData(obj, val)
      obj.DropDown.ItemsData = val;
    end
    
    function val = get.Items(obj)
      val = obj.DropDown.Items;
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.DropDown.ValueChangedFcn = val;
    end
  end
end