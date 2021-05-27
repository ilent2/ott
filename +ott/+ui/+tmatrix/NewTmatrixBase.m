classdef (Abstract) NewTmatrixBase < ott.ui.support.AppTopLevel ...
    & ott.ui.support.AppProducer ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Base class for beam creation application windows.
%
% This class is not intended to be instantiated directly.
% The T-matrix is stored internally and written to the matlab workspace
% if a variable name is given for the shape.  To access the internal
% instance use:
%
%   app = ott.ui.tmatrix.<name-of-your-app>()
%   tmatrix = app.Data

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    windowSize = [300, 420];
  end
  
  properties (Access=public)
    MainGrid             matlab.ui.container.GridLayout
    ExtraGrid            matlab.ui.container.GridLayout
    ShapeName            ott.ui.support.VariableDropDown
    WavelengthSpinner    ott.ui.support.LabeledSpinner
    RelativeIndexSpinner ott.ui.support.LabeledSpinner
  end
  
  properties (Access=protected)
    wwidth = 140
  end
  
  methods (Access=protected)
    function setDefaultValues(app)
      app.VariableName.Value = '';
      app.ShapeName.Value = '';
      app.WavelengthSpinner.Value = 1e-6;
      app.RelativeIndexSpinner.Value = 1.2;
      app.UpdateButton.Level = 0;
      app.UpdateButton.clearErrors();
    end
    
    function createMainComponents(app)
      
      % Create grid
      app.MainGrid = uigridlayout(app.UIFigure, [4, 1]);
      app.MainGrid.RowHeight = repmat({32}, 1, 6);
      app.MainGrid.RowHeight{end-1} = '1x';
      app.MainGrid.ColumnWidth = {'1x'};
      app.MainGrid.ColumnSpacing = 1;
      app.MainGrid.RowSpacing = 1;
      
      % Output variable entry
      app.VariableName = ott.ui.support.OutputVariableEntry(app.MainGrid, ...
        'wwidth', app.wwidth);
      app.VariableName.Layout.Row = 1;
      app.VariableName.Layout.Column = 1;
      
      % Shape input
      app.ShapeName = ott.ui.support.VariableDropDown(app.MainGrid, ...
        'wwidth', app.wwidth, 'filter', 'ott.shape.Shape');
      app.ShapeName.Layout.Row = 2;
      app.ShapeName.Layout.Column = 1;
      app.registerRefreshInput(app.ShapeName);
      
      % Wavelength spinner
      app.WavelengthSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'wwidth', app.wwidth, 'label', 'Wavelength');
      app.WavelengthSpinner.Layout.Row = 3;
      app.WavelengthSpinner.Layout.Column = 1;
      app.WavelengthSpinner.Step = 1e-6;
      app.WavelengthSpinner.Limits = [0,Inf];
      app.WavelengthSpinner.LowerLimitInclusive = false;
      app.WavelengthSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Index relative spinner
      app.RelativeIndexSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'wwidth', app.wwidth, 'label', 'Rel. Index');
      app.RelativeIndexSpinner.Layout.Row = 4;
      app.RelativeIndexSpinner.Layout.Column = 1;
      app.RelativeIndexSpinner.Step = 0.1;
      app.RelativeIndexSpinner.Limits = [0,Inf];
      app.RelativeIndexSpinner.LowerLimitInclusive = false;
      app.RelativeIndexSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Create grid
      app.ExtraGrid = uigridlayout(app.MainGrid);
      app.ExtraGrid.Layout.Row = 5;
      app.ExtraGrid.Layout.Column = 1;
      app.ExtraGrid.Padding = [0, 0, 0, 0];
      app.ExtraGrid.ColumnWidth = {'1x'};
      app.ExtraGrid.ColumnSpacing = 1;
      app.ExtraGrid.RowSpacing = 1;
      
      % Update button
      app.UpdateButton = ott.ui.support.UpdateWithProgress(app.MainGrid, ...
        'wwidth', app.wwidth);
      app.UpdateButton.Layout.Row = 6;
      app.UpdateButton.Layout.Column = 1;
      
    end
  end
  
  methods
    function app = NewTmatrixBase()
      
      % Call window constructor first to create widgets
      app = app@ott.ui.support.AppTopLevel();
      
      % Then call AppProducer to connect production widgets
      app = app@ott.ui.support.AppProducer();
      
    end
  end
end