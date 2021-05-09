classdef PmParaxial < ott.ui.beam.NewBeamBase ...
    & ott.ui.support.GenerateCodeMenu
% Generate a beam using paraxial point matching

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'PmParaxial';

    nameText = 'Paraxial Pointmatched Beam';

    aboutText = ['Generates a beam using paraxial point matching from' ...
      ' a given E field profile.'];
    
    helpText = {ott.ui.beam.PmParaxial.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
  end
  
  properties (Access=protected)
    VariableNameEditFieldLabel      matlab.ui.control.Label
    VariableNameEditField           matlab.ui.control.EditField
    GenerateButton                  matlab.ui.control.Button
    WavelengthvacuumEditFieldLabel  matlab.ui.control.Label
    WavelengthvacuumEditField       matlab.ui.control.EditField
    NmaxSpinnerLabel                matlab.ui.control.Label
    NmaxSpinner                     matlab.ui.control.Spinner
    AutoCheckBox                    matlab.ui.control.CheckBox
    IndexmediumEditFieldLabel       matlab.ui.control.Label
    IndexmediumEditField            matlab.ui.control.EditField
    NASpinnerLabel                  matlab.ui.control.Label
    NASpinner                       matlab.ui.control.Spinner
    MappingModeDropDownLabel        matlab.ui.control.Label
    MappingModeDropDown             matlab.ui.control.DropDown
    Gauge                           matlab.ui.control.LinearGauge
    EfarfieldmatrixLabel            matlab.ui.control.Label
    EfarfieldmatrixDropDown         matlab.ui.control.DropDown
    PolarisationEditFieldLabel      matlab.ui.control.Label
    PolarisationEditField           matlab.ui.control.EditField
    WarningNAlargerthanIndexmediumLabel  matlab.ui.control.Label
  end
  
  methods (Access=protected)
    
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function data = generateData(app)
      data = []; % TODO
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);

      % Create VariableNameEditFieldLabel
      app.VariableNameEditFieldLabel = uilabel(app.LeftPanel);
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

      % Create WavelengthvacuumEditFieldLabel
      app.WavelengthvacuumEditFieldLabel = uilabel(app.LeftPanel);
      app.WavelengthvacuumEditFieldLabel.Position = [14 259 122 22];
      app.WavelengthvacuumEditFieldLabel.Text = 'Wavelength (vacuum)';

      % Create WavelengthvacuumEditField
      app.WavelengthvacuumEditField = uieditfield(app.LeftPanel, 'text');
      app.WavelengthvacuumEditField.Position = [148 259 100 22];
      app.WavelengthvacuumEditField.Value = '1.0';

      % Create NmaxSpinnerLabel
      app.NmaxSpinnerLabel = uilabel(app.LeftPanel);
      app.NmaxSpinnerLabel.Position = [14 139 37 22];
      app.NmaxSpinnerLabel.Text = 'Nmax';

      % Create NmaxSpinner
      app.NmaxSpinner = uispinner(app.LeftPanel);
      app.NmaxSpinner.Limits = [1 Inf];
      app.NmaxSpinner.RoundFractionalValues = 'on';
      app.NmaxSpinner.Enable = 'off';
      app.NmaxSpinner.Position = [147 139 100 22];
      app.NmaxSpinner.Value = 30;

      % Create AutoCheckBox
      app.AutoCheckBox = uicheckbox(app.LeftPanel);
      app.AutoCheckBox.ValueChangedFcn = createCallbackFcn(app, @AutoCheckBoxValueChanged, true);
      app.AutoCheckBox.Text = 'Auto';
      app.AutoCheckBox.Position = [76 139 47 22];
      app.AutoCheckBox.Value = true;

      % Create IndexmediumEditFieldLabel
      app.IndexmediumEditFieldLabel = uilabel(app.LeftPanel);
      app.IndexmediumEditFieldLabel.Position = [14 222 89 22];
      app.IndexmediumEditFieldLabel.Text = 'Index (medium)';

      % Create IndexmediumEditField
      app.IndexmediumEditField = uieditfield(app.LeftPanel, 'text');
      app.IndexmediumEditField.ValueChangedFcn = createCallbackFcn(app, @IndexmediumEditFieldValueChanged, true);
      app.IndexmediumEditField.Position = [148 222 100 22];
      app.IndexmediumEditField.Value = '1.0';

      % Create NASpinnerLabel
      app.NASpinnerLabel = uilabel(app.LeftPanel);
      app.NASpinnerLabel.Position = [14 185 25 22];
      app.NASpinnerLabel.Text = 'NA';

      % Create NASpinner
      app.NASpinner = uispinner(app.LeftPanel);
      app.NASpinner.Step = 0.1;
      app.NASpinner.LowerLimitInclusive = 'off';
      app.NASpinner.Limits = [0 Inf];
      app.NASpinner.ValueChangedFcn = createCallbackFcn(app, @IndexmediumEditFieldValueChanged, true);
      app.NASpinner.Position = [148 185 100 22];
      app.NASpinner.Value = 0.9;

      % Create MappingModeDropDownLabel
      app.MappingModeDropDownLabel = uilabel(app.LeftPanel);
      app.MappingModeDropDownLabel.Position = [14 343 85 22];
      app.MappingModeDropDownLabel.Text = 'Mapping Mode';

      % Create MappingModeDropDown
      app.MappingModeDropDown = uidropdown(app.LeftPanel);
      app.MappingModeDropDown.Items = {'sin(theta)', 'tan(theta)', 'theta'};
      app.MappingModeDropDown.ItemsData = {'sintheta', 'tantheta', 'theta'};
      app.MappingModeDropDown.Position = [147 343 100 22];
      app.MappingModeDropDown.Value = 'sintheta';

      % Create Gauge
      app.Gauge = uigauge(app.LeftPanel, 'linear');
      app.Gauge.Enable = 'off';
      app.Gauge.FontSize = 8;
      app.Gauge.Position = [0 10 133 29];

      % Create EfarfieldmatrixLabel
      app.EfarfieldmatrixLabel = uilabel(app.LeftPanel);
      app.EfarfieldmatrixLabel.Position = [14 375 101 22];
      app.EfarfieldmatrixLabel.Text = 'E far-field (matrix)';

      % Create EfarfieldmatrixDropDown
      app.EfarfieldmatrixDropDown = uidropdown(app.LeftPanel);
      app.EfarfieldmatrixDropDown.Items = {};
      app.EfarfieldmatrixDropDown.Editable = 'on';
      app.EfarfieldmatrixDropDown.BackgroundColor = [1 1 1];
      app.EfarfieldmatrixDropDown.Position = [147 375 100 22];
      app.EfarfieldmatrixDropDown.Value = {};

      % Create PolarisationEditFieldLabel
      app.PolarisationEditFieldLabel = uilabel(app.LeftPanel);
      app.PolarisationEditFieldLabel.Position = [14 306 68 22];
      app.PolarisationEditFieldLabel.Text = 'Polarisation';

      % Create PolarisationEditField
      app.PolarisationEditField = uieditfield(app.LeftPanel, 'text');
      app.PolarisationEditField.Position = [147 306 100 22];
      app.PolarisationEditField.Value = '[1, 1i]';

      % Create WarningNAlargerthanIndexmediumLabel
      app.WarningNAlargerthanIndexmediumLabel = uilabel(app.LeftPanel);
      app.WarningNAlargerthanIndexmediumLabel.Position = [14 60 221 22];
      app.WarningNAlargerthanIndexmediumLabel.Text = 'Warning: NA larger than Index (medium)';
      
    end
  end
  
  methods (Access=public)
    function app=PmParaxial()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.NewBeamBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end