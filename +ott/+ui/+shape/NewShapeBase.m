classdef (Abstract) NewShapeBase < ott.ui.support.AppTwoColumn
% Base class for shape creation application windows.
%
% This class is not intended to be instantiated directly.
% The shape is stored internally and/or written to the matlab workspace
% if a variable name is given for the shape.  To access the internal
% instance use:
%
%   app = ott.ui.shape.<name-of-your-app>()
%   shape = app.shape

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Maybe merge this with NewBeamBase?

  properties (Constant)
    windowSize = [640, 420];
  end

  properties (SetAccess=protected)
    shape       % Internal representation of the shape
  end
  
  properties (Access=public)
    
    % Left panel
    MainGrid              matlab.ui.container.GridLayout
    ExtraGrid             matlab.ui.container.GridLayout
    VariableName          ott.ui.support.OutputVariableEntry
    OffsetXyzSpinners     ott.ui.support.XyzSpinners
    RotationXyzSpinners   ott.ui.support.XyzSpinners
    ShowPreviewCheckBox   matlab.ui.control.CheckBox
    UpdateButton          ott.ui.support.UpdateCheckButton
    
    % Right pannel
    LoadingText           matlab.ui.control.Label
    PreviewUIAxes         matlab.ui.control.UIAxes
    
  end
  
  methods (Access=protected, Abstract)
    generateShape(app)
  end
  
  methods (Access=protected)
    
    function startupFcn(app)
      app.updateShapePreview();
    end
    
    function updateShapePreview(app)
      
      if isempty(app.shape)
        return;
      end
      
      app.LoadingText.Visible = 'on';
      drawnow nocallbacks;
      
      % Update shape preview
      app.shape.surf('axes', app.PreviewUIAxes);
      
      app.LoadingText.Visible = 'off';
      
    end
    
    function setDefaultValues(app, ~)
      % Set default values to window fields
      
      app.VariableName.Value = '';
      app.OffsetXyzSpinners.Value = [0, 0, 0];
      app.RotationXyzSpinners.Value = [0, 0, 0];
      app.ShowPreviewCheckBox.Value = true;
      app.UpdateButton.Value = true;
      
      app.valueChangedCb();
      
    end
    
    function valueChangedCb(app, ~)
      % Regenerate data if auto-update is enabled
      if app.UpdateButton.Value
        app.updateCb();
      end
    end
    
    function updateCb(app, ~)
      % Called when a value is changed or when update is clicked
      
      % Generate new beam
      app.shape = app.generateShape();
      
      % Write to workspace (if requested, not needed if only preview)
      if ~isempty(app.VariableName.Value)
        app.VariableName.WriteVariable(app.shape);
      end
      
      % Generate beam preview
      app.showPreviewChangedCb();
    end
    
    function showPreviewChangedCb(app, ~)
      % Generate a new preview (without generating new beam)
      if app.ShowPreviewCheckBox.Value
        app.updateShapePreview();
      end
    end
    
    function nameChangedCb(app, ~)
      % Change output variable, dont regenerate data
      
      % Write to workspace (if requested, not needed if only preview)
      if ~isempty(app.VariableName.Value)
        app.VariableName.WriteVariable(app.beam);
      end
    end
    
    function createLeftComponents(app)
      
      % Create grid
      app.MainGrid = uigridlayout(app.LeftPanel);
      app.MainGrid.RowHeight = repmat({32}, 1, 6);
      app.MainGrid.RowHeight{end-2} = '1x';
      app.MainGrid.ColumnWidth = {230};
      app.MainGrid.ColumnSpacing = 1;
      app.MainGrid.RowSpacing = 1;
      
      % Variable name entry
      app.VariableName = ott.ui.support.OutputVariableEntry(app.MainGrid);
      app.VariableName.Layout.Row = 1;
      app.VariableName.Layout.Column = 1;
      app.VariableName.ValueChangedFcn = createCallbackFcn(app, ...
          @nameChangedCb, true);
        
      % Offset entry
      app.OffsetXyzSpinners = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Offset');
      app.OffsetXyzSpinners.Layout.Row = 2;
      app.OffsetXyzSpinners.Layout.Column = 1;
      app.OffsetXyzSpinners.Step = 1e-7;
      app.OffsetXyzSpinners.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Rotation entry
      app.RotationXyzSpinners = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Rotation');
      app.RotationXyzSpinners.Layout.Row = 3;
      app.RotationXyzSpinners.Layout.Column = 1;
      app.RotationXyzSpinners.Step = 10;
      app.RotationXyzSpinners.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Create grid
      app.ExtraGrid = uigridlayout(app.MainGrid);
      app.ExtraGrid.Padding = [0, 0, 0, 0];
      app.ExtraGrid.ColumnWidth = {'1x'};
      app.ExtraGrid.ColumnSpacing = 1;
      app.ExtraGrid.RowSpacing = 1;
      app.ExtraGrid.Layout.Row = 4;
      app.ExtraGrid.Layout.Column = 1;
      
      % Preview checkbox
      app.ShowPreviewCheckBox = uicheckbox(app.MainGrid);
      app.ShowPreviewCheckBox.Text = 'Show Preview';
      app.ShowPreviewCheckBox.Layout.Row = 5;
      app.ShowPreviewCheckBox.Layout.Column = 1;
      app.ShowPreviewCheckBox.Value = true;
      app.ShowPreviewCheckBox.ValueChangedFcn = createCallbackFcn(app, ...
          @showPreviewChangedCb, true);
      
      % Update button
      app.UpdateButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.UpdateButton.Layout.Row = 6;
      app.UpdateButton.Layout.Column = 1;
      addlistener(app.UpdateButton, "UpdateCalled", @(h,e) app.updateCb(e));
        
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
      
    end
  end
end