classdef UpdateCheckButton < ott.ui.support.GridWidget ...
    & ott.ui.support.UpdateButtonBase
% A check box and button for enabling/disabling auto-update and update.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  events (NotifyAccess = protected)
    UpdateCalled    % Emitted when button clicked or auto-update enabled
  end

  properties (Dependent)
    AutoUpdate           % Current auto-update value
  end
  
  properties
    CheckBox    matlab.ui.control.CheckBox
    Button      matlab.ui.control.Button
  end
  
  methods (Access=private)
    function onCheckBoxChanged(obj, evt)
      
      % Enable/disable button
      if strcmpi(obj.Button.Enable, 'on')
        obj.Button.Enable = 'off';
      else
        obj.Button.Enable = 'on';
      end
      
      % Emit UpdateCalled event (lazy: reuse input event...)
      if obj.Value
        notify(obj, "UpdateCalled", evt);
      end
      
    end
    
    function ButtonPushedCb(obj, evt)
      
      % Emit UpdateCalled event (lazy: reuse input event...)
      notify(obj, "UpdateCalled", evt);
    end
  end
  
  methods
    function obj = UpdateCheckButton(parent, varargin)
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      % Configure grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {'1x', 70};
      obj.Grid.ColumnSpacing = 10;
      obj.Grid.RowSpacing = 1;
      
      p = inputParser;
      p.addParameter('position', [0, 0]);
      p.parse(varargin{:});
      
      % Create AutoupdateCheckBox
      obj.CheckBox = uicheckbox(obj.Grid);
      obj.CheckBox.Text = 'Auto-update';
      obj.CheckBox.Layout.Column = 1;
      obj.CheckBox.Layout.Row = 1;
      obj.CheckBox.Value = true;
      obj.CheckBox.ValueChangedFcn = @(h,e) obj.onCheckBoxChanged(e);
      
      % Create UpdateButton
      obj.Button = uibutton(obj.Grid, 'push');
      obj.Button.Enable = 'off';
      obj.Button.Layout.Column = 2;
      obj.Button.Layout.Row = 1;
      obj.Button.Text = 'Update';
      obj.Button.ButtonPushedFcn = @(h,e) obj.ButtonPushedCb(e);
      
    end
  end
  
  methods
    function val = get.AutoUpdate(obj)
      val = obj.CheckBox.Value;
    end
    
    function set.AutoUpdate(obj, val)
      obj.CheckBox.Value = val;
    end
  end
end