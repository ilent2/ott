classdef JonesPolarisationEntry < ott.ui.support.GridWidget
% Entry for a Jones polarisation vector.
%
% This is currently just a text entry field.  In future it might be
% nice to change this to some kind of pop-up window.  But we really
% don't have time to implement that now.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    Label         matlab.ui.control.Label
    EditField     matlab.ui.control.EditField
  end
  
  properties (Dependent)
    Value
    ValueChangedFcn
  end
  
  methods
    function obj = JonesPolarisationEntry(parent, varargin)
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      p = inputParser;
      p.parse(varargin{:});
      
      % Create grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {'1x', 100};
      obj.Grid.ColumnSpacing = 1;
      obj.Grid.RowSpacing = 1;
      
      % Create label
      obj.Label = uilabel(obj.Grid);
      obj.Label.HorizontalAlignment = 'left';
      obj.Label.Layout.Column = 1;
      obj.Label.Layout.Row = 1;
      obj.Label.Text = 'Polarisation';

      % Create ScatteredBeamEditField
      obj.EditField = uieditfield(obj.Grid, 'text');
      obj.EditField.Layout.Column = 2;
      obj.EditField.Layout.Row = 1;
      obj.EditField.Value = '[1, 1i]';
      
    end
  end
  
  methods
    function val = get.Value(obj)
      val = evalin('base', obj.EditField.Value);
    end

    function set.Value(obj, val)
      obj.EditField.Value = ['[', num2str(val) ']'];
    end
    
    function set.ValueChangedFcn(obj, val)
      obj.EditField.ValueChangedFcn = val;
    end
  end
end