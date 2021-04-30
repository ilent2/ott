classdef UpdateCheckButton < handle
% A check box and button for enabling/disabling auto-update and update.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  events (NotifyAccess = protected)
    UpdateCalled    % Emitted when button clicked or auto-update enabled
  end

  properties (Dependent)
    Value           % Current auto-update value
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
      if strcmpi(obj.Value, 'on')
        notify(obj, "UpdateCalled", evt);
      end
      
    end
  end
  
  methods
    function obj = UpdateCheckButton(panel, varargin)
      
      p = inputParser;
      p.addParameter('position', [0, 0]);
      p.parse(varargin{:});
      
      height = 22;
      
      % Create AutoupdateCheckBox
      obj.CheckBox = uicheckbox(panel);
      obj.CheckBox.Text = 'Auto-update';
      obj.CheckBox.Position = [p.Results.position 87 height];
      obj.CheckBox.Value = true;
      obj.CheckBox.ValueChangedFcn = @(h,e) obj.onCheckBoxChanged(e);
      
      % Create UpdateButton
      obj.Button = uibutton(panel, 'push');
      obj.Button.Enable = 'off';
      obj.Button.Position = [([150, 0] + p.Results.position), 84 height];
      obj.Button.Text = 'Update';
      
    end
  end
  
  methods
    function val = get.Value(obj)
      val = obj.CheckBox.Value;
    end
    
    function set.Value(obj, val)
      obj.CheckBox.Value = val;
    end
  end
  
end