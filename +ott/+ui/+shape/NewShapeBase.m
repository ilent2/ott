classdef (Abstract) NewShapeBase < ott.ui.support.AppTwoColumn ...
    & ott.ui.support.AppProducer
% Base class for shape creation application windows.
%
% This class is not intended to be instantiated directly.
% The shape is stored internally and/or written to the matlab workspace
% if a variable name is given for the shape.  To access the internal
% instance use:
%
%   app = ott.ui.shape.<name-of-your-app>()
%   shape = app.Data

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    windowSize = [640, 420];
  end
  
  properties (Access=public)
    
    % Left panel
    MainGrid              matlab.ui.container.GridLayout
    ExtraGrid             matlab.ui.container.GridLayout
    OffsetXyzSpinners     ott.ui.support.XyzSpinners
    RotationXyzSpinners   ott.ui.support.XyzSpinners
    LocalViewCheckBox     matlab.ui.control.CheckBox
    
    % Right pannel
    LoadingText           matlab.ui.control.Label
    PreviewUIAxes         matlab.ui.control.UIAxes
    
  end
  
  methods (Access=protected)
    function setDefaultValues(app, update)
      % Set default values to window fields
      
      app.VariableName.Value = '';
      app.OffsetXyzSpinners.Value = [0, 0, 0];
      app.RotationXyzSpinners.Value = [0, 0, 0];
      app.LocalViewCheckBox.Value = true;
      app.ShowPreviewCheckBox.Value = true;
      app.UpdateButton.AutoUpdate = true;
      
      % Only call update if requested
      if nargin == 1 || update == true
        app.update();
      end
    end
    
    function updatePreview(app)
      
      app.LoadingText.Visible = 'on';
      drawnow nocallbacks;
      
      if app.LocalViewCheckBox.Value
        origin = 'local';
      else
        origin = 'global';
      end
      
      % Update shape preview
      app.Data.surf('axes', app.PreviewUIAxes, 'origin', origin);
      
      app.LoadingText.Visible = 'off';
      
    end
    
    function createLeftComponents(app)
      
      % Create grid
      app.MainGrid = uigridlayout(app.LeftPanel);
      app.MainGrid.RowHeight = repmat({32}, 1, 7);
      app.MainGrid.RowHeight{end-3} = '1x';
      app.MainGrid.ColumnWidth = {230};
      app.MainGrid.ColumnSpacing = 1;
      app.MainGrid.RowSpacing = 1;
      
      % Variable name entry
      app.VariableName = ott.ui.support.OutputVariableEntry(app.MainGrid);
      app.VariableName.Layout.Row = 1;
      app.VariableName.Layout.Column = 1;
        
      % Offset entry
      app.OffsetXyzSpinners = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Offset');
      app.OffsetXyzSpinners.Layout.Row = 2;
      app.OffsetXyzSpinners.Layout.Column = 1;
      app.OffsetXyzSpinners.Step = 1e-7;
      app.OffsetXyzSpinners.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Rotation entry
      app.RotationXyzSpinners = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Rotation');
      app.RotationXyzSpinners.Layout.Row = 3;
      app.RotationXyzSpinners.Layout.Column = 1;
      app.RotationXyzSpinners.Step = 10;
      app.RotationXyzSpinners.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Create grid
      app.ExtraGrid = uigridlayout(app.MainGrid);
      app.ExtraGrid.Padding = [0, 0, 0, 0];
      app.ExtraGrid.ColumnWidth = {'1x'};
      app.ExtraGrid.ColumnSpacing = 1;
      app.ExtraGrid.RowSpacing = 1;
      app.ExtraGrid.Layout.Row = 4;
      app.ExtraGrid.Layout.Column = 1;
      
      % Local view checkbox
      app.LocalViewCheckBox = uicheckbox(app.MainGrid);
      app.LocalViewCheckBox.Text = 'Local View';
      app.LocalViewCheckBox.Layout.Row = 5;
      app.LocalViewCheckBox.Layout.Column = 1;
      app.LocalViewCheckBox.ValueChangedFcn = @(~,~) app.updatePreview();
      
      % Preview checkbox
      app.ShowPreviewCheckBox = uicheckbox(app.MainGrid);
      app.ShowPreviewCheckBox.Text = 'Show Preview';
      app.ShowPreviewCheckBox.Layout.Row = 6;
      app.ShowPreviewCheckBox.Layout.Column = 1;
      
      % Update button
      app.UpdateButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.UpdateButton.Layout.Row = 7;
      app.UpdateButton.Layout.Column = 1;
        
    end
    
    function createRightComponents(app)
      
      % Create UIAxes
      app.PreviewUIAxes = uiaxes(app.RightPanel);
      title(app.PreviewUIAxes, 'Preview')
      xlabel(app.PreviewUIAxes, '')
      ylabel(app.PreviewUIAxes, '')
      view(app.PreviewUIAxes, [-40, 30]);
      app.PreviewUIAxes.Position = [10 10 373 app.windowSize(2)-20];
      
      % Create loading text
      app.LoadingText = uilabel(app.RightPanel);
      app.LoadingText.Position = [150, app.windowSize(2)/2, 70, 22];
      app.LoadingText.Text = 'Loading...';
      app.LoadingText.BackgroundColor = app.UIFigure.Color;
      app.LoadingText.HorizontalAlignment = 'center';
      app.LoadingText.Visible = 'off';
      
    end
  end
  
  methods
    function app = NewShapeBase()
      
      % Call window constructor first to create widgets
      app = app@ott.ui.support.AppTwoColumn();
      
      % Then call AppProducer to connect production widgets
      app = app@ott.ui.support.AppProducer();
      
    end
  end
end