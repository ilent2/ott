classdef Gaussian < ott.ui.beam.AppBase
% Generate a simple beam representation and visualise.
%
% Supported beams:
%   - Gaussian
%   - Laguerre-Gaussian
%   - Hermite-Gaussian
%   - Plane Wave
%   - Mathieu
%   - Webber
%   - Bessel
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

% TODO: Should we split this interface into Gaussian, PlaneWave
% and Annular?

  properties (Constant)
    cnameText = 'Gaussian';

    nameText = 'Create a Gaussian Beam';

    aboutText = ['Generate a Gaussian or related beam.  Can be used to' ...
      ' create Gaussian, Laguerre-Gaussian, Hermite-Gaussian ' ...
      'and Ince-Gaussian beams.'];
    
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
      
      % Create VariableNameEditFieldLabel
      app.VariableNameEditFieldLabel = uilabel(app.LeftPanel);
      app.VariableNameEditFieldLabel.HorizontalAlignment = 'right';
      app.VariableNameEditFieldLabel.Position = [14 427 84 22];
      app.VariableNameEditFieldLabel.Text = 'Variable Name';

      % Create VariableNameEditField
      app.VariableNameEditField = uieditfield(app.LeftPanel, 'text');
      app.VariableNameEditField.Position = [148 427 100 22];
      app.VariableNameEditField.Value = 'Beam';

      % Create GenerateButton
      app.GenerateButton = uibutton(app.LeftPanel, 'push');
      app.GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateButtonPushed, true);
      app.GenerateButton.Position = [148 13 100 22];
      app.GenerateButton.Text = 'Generate';

      % Create RadialModeSpinnerLabel
      app.RadialModeSpinnerLabel = uilabel(app.LeftPanel);
      app.RadialModeSpinnerLabel.HorizontalAlignment = 'right';
      app.RadialModeSpinnerLabel.Position = [14 380 73 22];
      app.RadialModeSpinnerLabel.Text = 'Radial Mode';

      % Create RadialModeSpinner
      app.RadialModeSpinner = uispinner(app.LeftPanel);
      app.RadialModeSpinner.Limits = [0 Inf];
      app.RadialModeSpinner.RoundFractionalValues = 'on';
      app.RadialModeSpinner.Position = [148 380 100 22];

      % Create AzimuthalModeSpinnerLabel
      app.AzimuthalModeSpinnerLabel = uilabel(app.LeftPanel);
      app.AzimuthalModeSpinnerLabel.HorizontalAlignment = 'right';
      app.AzimuthalModeSpinnerLabel.Position = [14 343 91 22];
      app.AzimuthalModeSpinnerLabel.Text = 'Azimuthal Mode';

      % Create AzimuthalModeSpinner
      app.AzimuthalModeSpinner = uispinner(app.LeftPanel);
      app.AzimuthalModeSpinner.RoundFractionalValues = 'on';
      app.AzimuthalModeSpinner.Position = [148 343 100 22];

      % Create WavelengthvacuumEditFieldLabel
      app.WavelengthvacuumEditFieldLabel = uilabel(app.LeftPanel);
      app.WavelengthvacuumEditFieldLabel.HorizontalAlignment = 'right';
      app.WavelengthvacuumEditFieldLabel.Position = [14 225 122 22];
      app.WavelengthvacuumEditFieldLabel.Text = 'Wavelength (vacuum)';

      % Create WavelengthvacuumEditField
      app.WavelengthvacuumEditField = uieditfield(app.LeftPanel, 'text');
      app.WavelengthvacuumEditField.Position = [148 226 100 22];
      app.WavelengthvacuumEditField.Value = '1.0';

      % Create NmaxSpinnerLabel
      app.NmaxSpinnerLabel = uilabel(app.LeftPanel);
      app.NmaxSpinnerLabel.HorizontalAlignment = 'right';
      app.NmaxSpinnerLabel.Position = [14 114 37 22];
      app.NmaxSpinnerLabel.Text = 'Nmax';

      % Create NmaxSpinner
      app.NmaxSpinner = uispinner(app.LeftPanel);
      app.NmaxSpinner.Limits = [1 Inf];
      app.NmaxSpinner.RoundFractionalValues = 'on';
      app.NmaxSpinner.Enable = 'off';
      app.NmaxSpinner.Position = [148 115 100 22];
      app.NmaxSpinner.Value = 20;

      % Create AutoCheckBox
      app.AutoCheckBox = uicheckbox(app.LeftPanel);
      app.AutoCheckBox.ValueChangedFcn = createCallbackFcn(app, @AutoCheckBoxValueChanged, true);
      app.AutoCheckBox.Text = 'Auto';
      app.AutoCheckBox.Position = [89 114 47 22];
      app.AutoCheckBox.Value = true;

      % Create IndexmediumEditFieldLabel
      app.IndexmediumEditFieldLabel = uilabel(app.LeftPanel);
      app.IndexmediumEditFieldLabel.HorizontalAlignment = 'right';
      app.IndexmediumEditFieldLabel.Position = [14 188 89 22];
      app.IndexmediumEditFieldLabel.Text = 'Index (medium)';

      % Create IndexmediumEditField
      app.IndexmediumEditField = uieditfield(app.LeftPanel, 'text');
      app.IndexmediumEditField.ValueChangedFcn = createCallbackFcn(app, @IndexmediumEditFieldValueChanged, true);
      app.IndexmediumEditField.Position = [148 189 100 22];
      app.IndexmediumEditField.Value = '1.0';

      % Create NASpinnerLabel
      app.NASpinnerLabel = uilabel(app.LeftPanel);
      app.NASpinnerLabel.HorizontalAlignment = 'right';
      app.NASpinnerLabel.Position = [14 151 25 22];
      app.NASpinnerLabel.Text = 'NA';

      % Create NASpinner
      app.NASpinner = uispinner(app.LeftPanel);
      app.NASpinner.Step = 0.1;
      app.NASpinner.LowerLimitInclusive = 'off';
      app.NASpinner.Limits = [0 Inf];
      app.NASpinner.ValueChangedFcn = createCallbackFcn(app, @IndexmediumEditFieldValueChanged, true);
      app.NASpinner.Position = [148 152 100 22];
      app.NASpinner.Value = 0.9;

      % Create Gauge
      app.Gauge = uigauge(app.LeftPanel, 'linear');
      app.Gauge.Enable = 'off';
      app.Gauge.FontSize = 8;
      app.Gauge.Position = [1 10 133 29];

      % Create WarningNAlargerthanIndexmediumLabel
      app.WarningNAlargerthanIndexmediumLabel = uilabel(app.LeftPanel);
      app.WarningNAlargerthanIndexmediumLabel.Position = [10 61 221 22];
      app.WarningNAlargerthanIndexmediumLabel.Text = 'Warning: NA larger than Index (medium)';

      % Create PolarisationEditFieldLabel
      app.PolarisationEditFieldLabel = uilabel(app.LeftPanel);
      app.PolarisationEditFieldLabel.HorizontalAlignment = 'right';
      app.PolarisationEditFieldLabel.Position = [14 300 68 22];
      app.PolarisationEditFieldLabel.Text = 'Polarisation';

      % Create PolarisationEditField
      app.PolarisationEditField = uieditfield(app.LeftPanel, 'text');
      app.PolarisationEditField.Position = [148 301 100 22];
      app.PolarisationEditField.Value = '[1, 1i]';

      % Create PowerEditFieldLabel
      app.PowerEditFieldLabel = uilabel(app.LeftPanel);
      app.PowerEditFieldLabel.HorizontalAlignment = 'right';
      app.PowerEditFieldLabel.Position = [14 267 39 22];
      app.PowerEditFieldLabel.Text = 'Power';

      % Create PowerEditField
      app.PowerEditField = uieditfield(app.LeftPanel, 'text');
      app.PowerEditField.Position = [148 268 100 22];
      app.PowerEditField.Value = '1.0';
      
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