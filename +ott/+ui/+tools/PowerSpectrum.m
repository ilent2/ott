classdef PowerSpectrum < ott.ui.support.AppTwoColumn ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Generate plot of power spectrum from position time series data
%
% This app cab be launched from the app launcher via Tools -> PowerSpectrum
% or instantiated on the command line with:
%
%   ott.ui.tools.PowerSpectrum()

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Add support for input arguments to class

  properties (Constant)
    cnameText = 'PowerSpectrum';

    nameText = 'Generate Power Spectrum';

    aboutText = ['Generate plots of power spectrum from position time' ...
      'series data.'];
    
    helpText = {ott.ui.tools.PowerSpectrum.aboutText, ...
      ''};
    
    windowName = ott.ui.tools.PowerSpectrum.nameText;
    windowSize = [640, 400];
  end
  
  properties
    
    % Left panel
    MainGrid          matlab.ui.container.GridLayout
    XInputDropDown    ott.ui.support.VariableDropDown
    YInputDropDown    ott.ui.support.VariableDropDown
    XLabelEntry       ott.ui.support.LabeledTextEntry
    YLabelEntry       ott.ui.support.LabeledTextEntry
    UpdateButton      ott.ui.support.UpdateCheckButton
    
    % Right panel
    UIAxes            matlab.ui.control.UIAxes
  end
  
  methods (Access=protected)
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function setDefaultValues(app, ~)
      app.XInputDropDown.Value = '';
      app.YInputDropDown.Value = '';
      app.XLabelEntry.Value = 'Frequency [a.u.]';
      app.YLabelEntry.Value = 'Power [a.u.]';
      app.UpdateButton.AutoUpdate = true;
      
      % Update the graph
      app.UpdateCallback();
    end
    
    function ValueChangedCallback(app, ~)
      % Called when a widget changes
      if app.UpdateButton.AutoUpdate
        app.UpdateCallback();
      end
    end
    
    function UpdateCallback(app, ~)
      % Called when auto-update is enabled or the update button is clicked
     
      % TODO: Update graph content
      
      % Update axis labels
      xlabel(app.UIAxes, app.XLabelEntry.Value);
      ylabel(app.UIAxes, app.YLabelEntry.Value);
      
    end
    
    function createLeftComponents(app)
      
      % Main layout grid
      app.MainGrid = uigridlayout(app.LeftPanel);
      app.MainGrid.RowHeight = repmat({32}, 1, 6);
      app.MainGrid.RowHeight{end-1} = '1x';
      app.MainGrid.ColumnWidth = {'1x'};
      app.MainGrid.RowSpacing = 0;
      
      wwidth = 120;
      
      % T variable selection
      app.XInputDropDown = ott.ui.support.VariableDropDown(app.MainGrid, ...
        'label', 'Time data', 'wwidth', wwidth, 'filter', 'double');
      app.XInputDropDown.Layout.Row = 1;
      app.XInputDropDown.Layout.Column = 1;
      app.XInputDropDown.ValueChangedFcn = createCallbackFcn(app, ...
          @ValueChangedCallback, true);
      app.registerRefreshInput(app.XInputDropDown);
      
      % Y variable selection
      app.YInputDropDown = ott.ui.support.VariableDropDown(app.MainGrid, ...
        'label', 'Value data', 'wwidth', wwidth, 'filter', 'double');
      app.YInputDropDown.Layout.Row = 2;
      app.YInputDropDown.Layout.Column = 1;
      app.YInputDropDown.ValueChangedFcn = createCallbackFcn(app, ...
          @ValueChangedCallback, true);
      app.registerRefreshInput(app.YInputDropDown);
      
      app.XLabelEntry = ott.ui.support.LabeledTextEntry(app.MainGrid, ...
        'label', 'Horz. Label', 'wwidth', wwidth);
      app.XLabelEntry.Layout.Row = 3;
      app.XLabelEntry.Layout.Column = 1;
      app.XLabelEntry.ValueChangedFcn = createCallbackFcn(app, ...
          @ValueChangedCallback, true);
      
      app.YLabelEntry = ott.ui.support.LabeledTextEntry(app.MainGrid, ...
        'label', 'Vert. Label', 'wwidth', wwidth);
      app.YLabelEntry.Layout.Row = 4;
      app.YLabelEntry.Layout.Column = 1;
      app.YLabelEntry.ValueChangedFcn = createCallbackFcn(app, ...
          @ValueChangedCallback, true);

      % Create CalculateButton
      app.UpdateButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.UpdateButton.Layout.Row = 6;
      app.UpdateButton.Layout.Column = 1;
      addlistener(app.UpdateButton, "UpdateCalled", ...
          @(h,e) app.UpdateCallback(e));
      
    end
    
    function createRightComponents(app)
      
      app.UIAxes = uiaxes(app.RightPanel);
      app.UIAxes.Position = [10, 10, 370, 380];
      
    end
  end
  
  methods (Access=public)
    function app=PowerSpectrum()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTwoColumn();
      
      if nargout == 0
        clear app;
      end
    end
  end
end