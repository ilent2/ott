classdef Gaussian < ott.ui.beam.NewBeamBase
% Generate a simple beam representation and visualise.
%
% Supported beams:
%   - Gaussian
%   - Laguerre-Gaussian
%   - Hermite-Gaussian
%
% Some of these beams might move to their own interface in a future
% release.
%
% This GUI can be launched from the launcher Beam -> Simple or by running:
%
%   ott.ui.beam.Simple

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Do we want a NA warning?

  properties (Constant)
    cnameText = 'Gaussian';

    nameText = 'Create a Gaussian Beam';

    aboutText = ['Generate a Gaussian or related beam.  Can be used to' ...
      ' create Gaussian, Laguerre-Gaussian, Hermite-Gaussian.'];
    
    helpText = {ott.ui.beam.Gaussian.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.Gaussian.nameText;
    windowSize = [640, 420];
  end
  
  properties (Access=public)
    VariableNameEditFieldLabel      matlab.ui.control.Label
    VariableNameEditField           matlab.ui.control.EditField
    GenerateButton                  matlab.ui.control.Button
    RadialModeSpinnerLabel          matlab.ui.control.Label
    RadialModeSpinner               matlab.ui.control.Spinner
    AzimuthalModeSpinnerLabel       matlab.ui.control.Label
    AzimuthalModeSpinner            matlab.ui.control.Spinner
    WavelengthvacuumEditFieldLabel  matlab.ui.control.Label
    WavelengthvacuumEditField       matlab.ui.control.EditField
    NmaxSpinnerLabel                matlab.ui.control.Label
    NmaxSpinner                     matlab.ui.control.Spinner
    AutoCheckBox                    matlab.ui.control.CheckBox
    IndexmediumEditFieldLabel       matlab.ui.control.Label
    IndexmediumEditField            matlab.ui.control.EditField
    NASpinnerLabel                  matlab.ui.control.Label
    NASpinner                       matlab.ui.control.Spinner
    Gauge                           matlab.ui.control.LinearGauge
    WarningNAlargerthanIndexmediumLabel  matlab.ui.control.Label
    PolarisationEditFieldLabel      matlab.ui.control.Label
    PolarisationEditField           matlab.ui.control.EditField
    PowerEditFieldLabel             matlab.ui.control.Label
    PowerEditField                  matlab.ui.control.EditField
  end
  
  methods (Access=protected)
    function startupFcn(app)
    end
    
    function createRightComponents(app)
      app.createBeamPreview();
    end
    
    function createLeftComponents(app)
      
      % Create grid
      app.MainGrid = uigridlayout(app.UIFigure);
      app.MainGrid.RowHeight = repmat({32}, 1, 15);
      app.MainGrid.RowHeight(end-2) = '1x';
      app.MainGrid.ColumnWidth = {230};
      app.MainGrid.ColumnSpacing = 1;
      app.MainGrid.RowSpacing = 1;
      
      % Variable name entry
      app.VariableName = ott.ui.support.OutputVariableEntry(app.MainGrid);
      app.VariableName.Layout.Row = 1;
      app.VariableName.Layout.Column = 1;
      app.VariableName.ValueChangedFcn = createCallbackFcn(app, ...
          @nameChangedCb, true);
        
      % Beam type dropdown
      app.BeamTypeDropdown = ott.ui.support.LabeledDropdown(app.MainGrid);
      app.BeamTypeDropdown.Layout.Row = 2;
      app.BeamTypeDropdown.Layout.Column = 1;
      app.BeamTypeDropdown.ValueChangedFcn = createCallbackFcn(app, ...
          @typeChangedCb, true);
        
      % Lmode spinner
      app.LmodeSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'L mode');
      app.LmodeSpinner.Layout.Row = 3;
      app.LmodeSpinner.Layout.Column = 1;
      app.LmodeSpinner.Step = 1;
      app.LmodeSpinner.Limits = [-Inf,Inf];
      app.LmodeSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Pmode spinner
      app.PmodeSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'P mode');
      app.PmodeSpinner.Layout.Row = 4;
      app.PmodeSpinner.Layout.Column = 1;
      app.PmodeSpinner.Step = 1;
      app.PmodeSpinner.Limits = [0,Inf];
      app.PmodeSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Mmode spinner
      app.MmodeSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'M mode');
      app.MmodeSpinner.Layout.Row = 5;
      app.MmodeSpinner.Layout.Column = 1;
      app.MmodeSpinner.Step = 1;
      app.MmodeSpinner.Limits = [0,Inf];
      app.MmodeSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Nmode spinner
      app.NmodeSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'N mode');
      app.NmodeSpinner.Layout.Row = 6;
      app.NmodeSpinner.Layout.Column = 1;
      app.NmodeSpinner.Step = 1;
      app.NmodeSpinner.Limits = [0,Inf];
      app.NmodeSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Wavelength spinner
      app.WavelengthSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Wavelength');
      app.WavelengthSpinner.Layout.Row = 7;
      app.WavelengthSpinner.Layout.Column = 1;
      app.WavelengthSpinner.Step = 1e-7;
      app.WavelengthSpinner.Limits = [1e-9, Inf];
      app.WavelengthSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Refractive index spinner
      app.IndexSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Refractive Index');
      app.IndexSpinner.Layout.Row = 8;
      app.IndexSpinner.Layout.Column = 1;
      app.IndexSpinner.Step = 0.1;
      app.IndexSpinner.Limits = [0.1, Inf];
      app.IndexSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
        
      % NA spinner
      app.NaSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Numerical Aperture');
      app.NaSpinner.Layout.Row = 9;
      app.NaSpinner.Layout.Column = 1;
      app.NaSpinner.Step = 0.1;
      app.NaSpinner.LowerLimitInclusive = 'off';
      app.NaSpinner.Limits = [0.0, Inf];
      app.NaSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
        
      % Polarisation jones vector entry
      app.PolarisationEntry = ott.ui.support.JonesPolarisationEntry(...
          app.MainGrid);
      app.PolarisationEntry.Layout.Row = 10;
      app.PolarisationEntry.Layout.Column = 1;
      app.PolarisationEntry.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Power spinner  
      app.PowerSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Power');
      app.PowerSpinner.Layout.Row = 11;
      app.PowerSpinner.Layout.Column = 1;
      app.PowerSpinner.Step = 0.1;
      app.PowerSpinner.LowerLimitInclusive = 'off';
      app.PowerSpinner.Limits = [0.0, Inf];
      app.PowerSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Direction entry
      app.RotationXyzSpinner = ott.ui.support.XyzSpinners(...
          app.MainGrid, 'label', 'Rotation');
      app.RotationXyzSpinner.Layout.Row = 12;
      app.RotationXyzSpinner.Layout.Column = 1;
      app.RotationXyzSpinner.Step = 10;
      app.RotationXyzSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
      % Preview checkbox
      app.ShowPreviewCheckBox = uicheckbox(app.MainGrid);
      app.ShowPreviewCheckBox.Text = 'Show Preview';
      app.ShowPreviewCheckBox.Layout.Row = 14;
      app.ShowPreviewCheckBox.Layout.Column = 1;
      app.ShowPreviewCheckBox.Value = true;
      app.ShowPreviewCheckBox.ValueChangedFcn = createCallbackFcn(app, ...
          @showPreviewChangedCb, true);
      
      % Update button
      app.UpdateButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.UpdateButton.Layout.Row = 15;
      app.UpdateButton.Layout.Column = 1;
      addlistener(app.UpdateButton, "UpdateCalled", @(h,e) app.updateCb(e));

    end
  end
  
  methods (Access=public)
    function app=Gaussian()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.AppBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end