classdef UpdateWithProgress < ott.ui.support.GridWidget ...
    & ott.ui.support.UpdateButtonBase
% Create an update button with an adjacent progress bar
  
% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    Gauge         matlab.ui.control.LinearGauge
    Button      matlab.ui.control.Button
  end
  
  properties (Dependent)
    AutoUpdate
    Level
  end
  
  methods (Access = private)
    function buttonPushedCb(obj)
      % Callback for button: dispatches notification
      
      % Emit UpdateCalled event
      notify(obj, "UpdateCalled");
    end
  end
  
  methods
    function obj = UpdateWithProgress(parent, varargin)
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      % Configure grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {'1x', 120};
      obj.Grid.ColumnSpacing = 10;
      obj.Grid.RowSpacing = 1;
      
      % Create Button
      obj.Button = uibutton(obj.Grid, 'push');
      obj.Button.Text = 'Update';
      obj.Button.Layout.Column = 1;
      obj.Button.Layout.Row = 1;
      obj.Button.ButtonPushedFcn = @(~,~) obj.buttonPushedCb();

      % Create Guage
      obj.Gauge = uigauge(obj.Grid, 'linear');
      obj.Gauge.Enable = 'off';
      obj.Gauge.FontSize = 8;
      obj.Gauge.Layout.Column = 2;
      obj.Gauge.Layout.Row = 1;
      
    end
    
    function clearErrors(obj)
      obj.Gauge.BackgroundColor = 'white';
    end
    
    function setError(obj)
      obj.Gauge.BackgroundColor = 'red';
    end
    
    function setWarning(obj)
      obj.Gauge.BackgroundColor = 'yellow';
    end
  end
  
  methods
    function val = get.AutoUpdate(~)
      val = false;
    end
    
    function set.AutoUpdate(~, ~)
      error('Can not set AutoUpdate for UpdateWithProgress widget');
    end
    
    function set.Level(app, val)
      app.Gauge.Value = val;
    end
  end
end