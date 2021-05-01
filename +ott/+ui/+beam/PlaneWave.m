classdef PlaneWave < ott.ui.beam.AppBase
% Generate a plane wave beam.
%
% This GUI can be launched from the launcher Beam -> PlaneWave or with:
%
%   ott.ui.beam.PlaneWave

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'PlaneWave';

    nameText = 'Create Plane Wave Beam';

    aboutText = ['Generate a plane wave beam.'];
    
    helpText = {ott.ui.beam.PlaneWave.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PlaneWave.nameText;
    windowSize = [640, 420];
  end
  
  properties (Access=public)
    MainGrid            matlab.ui.container.GridLayout
    VariableName        ott.ui.support.OutputVariableEntry
    WavelengthSpinner   ott.ui.support.LabeledSpinner
    IndexSpinner        ott.ui.support.LabeledSpinner
    PolarisationEntry   ott.ui.support.JonesPolarisationEntry
    RotationXyzSpinner  ott.ui.support.XyzSpinners
    ShowPreviewCheckBox matlab.ui.control.CheckBox
    UpdateButton        ott.ui.support.UpdateCheckButton
  end
  
  methods (Access=protected)
    function startupFcn(app)
      app.updateBeamPreview();
    end
    
    function setDefaultValues(app, evt)
      % Set default values to window fields
      
      app.VariableName.Value = '';
      app.WavelengthSpinner.Value = 1.0e-6;
      app.IndexSpinner.Value = 1.0;
      app.PolarisationEntry.Value = [1, 1i];
      app.RotationXyzSpinner.Value = [0, 0, 0];
      app.ShowPreviewCheckBox.Value = true;
      app.UpdateButton.Value = true;
      
      app.valueChangedCb();
      
    end
    
    function valueChangedCb(app, evt)
      % Regenerate data if auto-update is enabled
      
      if ~app.UpdateButton.Value
        return;
      end
      
      % Update the beam
      app.updateCb();
    end
    
    function updateCb(app, evt)
      
      rot3 = app.RotationXyzSpinner.Value;
      rotation = ott.utils.rotx(rot3(1))*ott.utils.roty(rot3(2))*ott.utils.rotz(rot3(3));
      
      % Generate new beam
      app.beam = ott.beam.PlaneWave(...
          'polarisation', app.PolarisationEntry.Value, ...
          'index_medium', app.IndexSpinner.Value, ...
          'wavelength0', app.WavelengthSpinner.Value, ...
          'rotation', rotation);
      
      % Generate beam preview
      app.showPreviewChangedCb();
    end
    
    function showPreviewChangedCb(app, evt)
      % Generate a new preview (without generating new beam)
      
      if app.ShowPreviewCheckBox.Value
        app.updateBeamPreview();
      end
    end
    
    function nameChangedCb(app, evt)
      % TODO: Change output variable, dont regenerate data
      warning('not implemented');
    end
    
    function createRightComponents(app)
      % This feels like kludge, might change
      app.createBeamPreview();
    end
    
    function createLeftComponents(app)
      
      % Create grid
      app.MainGrid = uigridlayout(app.UIFigure, [4, 1]);
      app.MainGrid.RowHeight = {32, 32, 32, 32, 32, '1x', 32, 32};
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
        
      % Polarisation jones vector entry
      app.PolarisationEntry = ott.ui.support.JonesPolarisationEntry(...
          app.MainGrid);
      app.PolarisationEntry.Layout.Row = 4;
      app.PolarisationEntry.Layout.Column = 1;
      app.PolarisationEntry.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Direction entry
      app.RotationXyzSpinner = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Rotation');
      app.RotationXyzSpinner.Layout.Row = 5;
      app.RotationXyzSpinner.Layout.Column = 1;
      app.RotationXyzSpinner.Step = 10;
      app.RotationXyzSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
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
  end
  
  methods (Access=public)
    function app=PlaneWave()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.AppBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end