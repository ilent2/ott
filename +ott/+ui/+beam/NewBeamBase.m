classdef (Abstract) NewBeamBase < ott.ui.support.AppTwoColumn

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    windowSize = [640, 420];
  end

  properties (SetAccess=protected)
    beam
  end
  
  properties (Access=public)
    
    % Left panel
    MainGrid              matlab.ui.container.GridLayout
    ExtraGrid             matlab.ui.container.GridLayout
    VariableName          ott.ui.support.OutputVariableEntry
    WavelengthSpinner     ott.ui.support.LabeledSpinner
    IndexSpinner          ott.ui.support.LabeledSpinner
    RotationXyzSpinner    ott.ui.support.XyzSpinners
    TranslationXyzSpinner ott.ui.support.XyzSpinners
    ShowPreviewCheckBox   matlab.ui.control.CheckBox
    UpdateButton          ott.ui.support.UpdateCheckButton
    
    % Right pannel
    LoadingText           matlab.ui.control.Label
    PreviewUIAxes         matlab.ui.control.UIAxes
  end
  
  methods (Access=protected, Abstract)
    generateBeam(app)
    generateCode(app)
  end
  
  methods (Access=protected)
    
    function startupFcn(app)
      
      app.updateBeamPreview();
    end
    
    function updateBeamPreview(app)
      
      if isempty(app.beam)
        return;
      end
      
      app.LoadingText.Visible = 'on';
      drawnow nocallbacks;
      
      app.beam.visNearfield('plot_axes', app.PreviewUIAxes, ...
        'axis', 'y', 'range', [1,1]*2e-6, 'field', 'Re(Ex)', ...
        'size', [60, 60]);
      
      app.LoadingText.Visible = 'off';
      
    end
    
    function setDefaultValues(app, evt)
      % Set default values to window fields
      
      app.VariableName.Value = '';
      app.WavelengthSpinner.Value = 1.0e-6;
      app.IndexSpinner.Value = 1.0;
      app.TranslationXyzSpinner.Value = [0, 0, 0];
      app.RotationXyzSpinner.Value = [0, 0, 0];
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
      app.beam = app.generateBeam();
      
      % Generate beam preview
      app.showPreviewChangedCb();
    end
    
    function showPreviewChangedCb(app, ~)
      % Generate a new preview (without generating new beam)
      if app.ShowPreviewCheckBox.Value
        app.updateBeamPreview();
      end
    end
    
    function nameChangedCb(app, ~)
      % TODO: Change output variable, dont regenerate data
      warning('not implemented');
    end
    
    function createLeftComponents(app)
      
      % Create grid
      app.MainGrid = uigridlayout(app.UIFigure);
      app.MainGrid.RowHeight = repmat({32}, 1, 8);
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
      
      % Wavelength spinner
      app.WavelengthSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Wavelength');
      app.WavelengthSpinner.Layout.Row = 2;
      app.WavelengthSpinner.Layout.Column = 1;
      app.WavelengthSpinner.Step = 1e-7;
      app.WavelengthSpinner.Limits = [1e-9, Inf];
      app.WavelengthSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Refractive index spinner
      app.IndexSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Refractive Index');
      app.IndexSpinner.Layout.Row = 3;
      app.IndexSpinner.Layout.Column = 1;
      app.IndexSpinner.Step = 0.1;
      app.IndexSpinner.Limits = [0.1, Inf];
      app.IndexSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Offset entry
      app.TranslationXyzSpinner = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Translation');
      app.TranslationXyzSpinner.Layout.Row = 4;
      app.TranslationXyzSpinner.Layout.Column = 1;
      app.TranslationXyzSpinner.Step = 1e-7;
      app.TranslationXyzSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Direction entry
      app.RotationXyzSpinner = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Rotation');
      app.RotationXyzSpinner.Layout.Row = 5;
      app.RotationXyzSpinner.Layout.Column = 1;
      app.RotationXyzSpinner.Step = 10;
      app.RotationXyzSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Create grid
      app.ExtraGrid = uigridlayout(app.MainGrid);
      app.ExtraGrid.Padding = [0, 0, 0, 0];
      app.ExtraGrid.ColumnWidth = {'1x'};
      app.ExtraGrid.ColumnSpacing = 1;
      app.ExtraGrid.RowSpacing = 1;
      app.ExtraGrid.Layout.Row = 6;
      app.ExtraGrid.Layout.Column = 1;
      
      % Preview checkbox
      app.ShowPreviewCheckBox = uicheckbox(app.MainGrid);
      app.ShowPreviewCheckBox.Text = 'Show Preview';
      app.ShowPreviewCheckBox.Layout.Row = 7;
      app.ShowPreviewCheckBox.Layout.Column = 1;
      app.ShowPreviewCheckBox.Value = true;
      app.ShowPreviewCheckBox.ValueChangedFcn = createCallbackFcn(app, ...
          @showPreviewChangedCb, true);
      
      % Update button
      app.UpdateButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.UpdateButton.Layout.Row = 8;
      app.UpdateButton.Layout.Column = 1;
      addlistener(app.UpdateButton, "UpdateCalled", @(h,e) app.updateCb(e));
      
    end
    
    function createRightComponents(app)
      
      % Create UIAxes
      app.PreviewUIAxes = uiaxes(app.RightPanel);
      title(app.PreviewUIAxes, 'Preview')
      xlabel(app.PreviewUIAxes, '')
      ylabel(app.PreviewUIAxes, '')
      app.PreviewUIAxes.XAxisLocation = 'origin';
      app.PreviewUIAxes.XTick = [];
      app.PreviewUIAxes.YAxisLocation = 'origin';
      app.PreviewUIAxes.YTick = [];
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