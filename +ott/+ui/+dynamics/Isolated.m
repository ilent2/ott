classdef Isolated < ott.ui.support.AppTwoColumn

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Isolated';

    nameText = 'Simulate Isolated Dynamics';

    aboutText = ['Simulate dynamics of an isolated particle.'];
    
    helpText = {ott.ui.dynamics.Isolated.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 420];
  end
  
  properties (Access=public)
    BeamDropDownLabel              matlab.ui.control.Label
    BeamDropDown                   matlab.ui.control.DropDown
    TmatrixDropDownLabel           matlab.ui.control.Label
    TmatrixDropDown                matlab.ui.control.DropDown
    ShapeDropDownLabel             matlab.ui.control.Label
    ShapeDropDown                  matlab.ui.control.DropDown
    SimulateButton                 matlab.ui.control.Button
    Gauge                          matlab.ui.control.LinearGauge
    OutputVariableEditFieldLabel   matlab.ui.control.Label
    OutputVariableEditField        matlab.ui.control.EditField
    TimestepSpinnerLabel           matlab.ui.control.Label
    TimestepSpinner                matlab.ui.control.Spinner
    DragEditFieldLabel             matlab.ui.control.Label
    DragEditField                  matlab.ui.control.EditField
    NumStepsSpinnerLabel           matlab.ui.control.Label
    NumStepsSpinner                matlab.ui.control.Spinner
    InitialPositionEditFieldLabel  matlab.ui.control.Label
    InitialPositionEditField       matlab.ui.control.EditField
    InitialRotationEditFieldLabel  matlab.ui.control.Label
    InitialRotationEditField       matlab.ui.control.EditField
    UIAxes                         matlab.ui.control.UIAxes
    UIAxes2                        matlab.ui.control.UIAxes
    DropDown                       matlab.ui.control.DropDown
  end
  
  methods (Access=protected)
    function startupFcn(app)
    end
    
    function createLeftComponents(app)
      
      lmargin = 10;

      % Create BeamDropDownLabel
      app.BeamDropDownLabel = uilabel(app.LeftPanel);
      app.BeamDropDownLabel.Position = [lmargin 390 37 22];
      app.BeamDropDownLabel.Text = 'Beam';

      % Create BeamDropDown
      app.BeamDropDown = uidropdown(app.LeftPanel);
      app.BeamDropDown.Items = {};
      app.BeamDropDown.Editable = 'on';
      app.BeamDropDown.BackgroundColor = [1 1 1];
      app.BeamDropDown.Position = [133 390 100 22];
      app.BeamDropDown.Value = {};

      % Create TmatrixDropDownLabel
      app.TmatrixDropDownLabel = uilabel(app.LeftPanel);
      app.TmatrixDropDownLabel.Position = [lmargin 360 49 22];
      app.TmatrixDropDownLabel.Text = 'Particle';

      % Create TmatrixDropDown
      app.TmatrixDropDown = uidropdown(app.LeftPanel);
      app.TmatrixDropDown.Items = {};
      app.TmatrixDropDown.Editable = 'on';
      app.TmatrixDropDown.BackgroundColor = [1 1 1];
      app.TmatrixDropDown.Position = [133 360 100 22];
      app.TmatrixDropDown.Value = {};

      % Create SimulateButton
      app.SimulateButton = uibutton(app.LeftPanel, 'push');
      app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateButtonPushed, true);
      app.SimulateButton.Position = [141 4 100 22];
      app.SimulateButton.Text = 'Simulate';

      % Create Gauge
      app.Gauge = uigauge(app.LeftPanel, 'linear');
      app.Gauge.Enable = 'off';
      app.Gauge.FontSize = 8;
      app.Gauge.Position = [1 1 133 29];

      % Create OutputVariableEditFieldLabel
      app.OutputVariableEditFieldLabel = uilabel(app.LeftPanel);
      app.OutputVariableEditFieldLabel.HorizontalAlignment = 'right';
      app.OutputVariableEditFieldLabel.Position = [lmargin 45 88 22];
      app.OutputVariableEditFieldLabel.Text = 'Output Variable';

      % Create OutputVariableEditField
      app.OutputVariableEditField = uieditfield(app.LeftPanel, 'text');
      app.OutputVariableEditField.ValueChangedFcn = createCallbackFcn(app, @OutputVariableEditFieldValueChanged, true);
      app.OutputVariableEditField.Position = [133 45 100 22];
      app.OutputVariableEditField.Value = 'SimulationData';

      % Create TimestepSpinnerLabel
      app.TimestepSpinnerLabel = uilabel(app.LeftPanel);
      app.TimestepSpinnerLabel.HorizontalAlignment = 'right';
      app.TimestepSpinnerLabel.Position = [lmargin 140 58 22];
      app.TimestepSpinnerLabel.Text = 'Time step';

      % Create TimestepSpinner
      app.TimestepSpinner = uispinner(app.LeftPanel);
      app.TimestepSpinner.Step = 1e-05;
      app.TimestepSpinner.LowerLimitInclusive = 'off';
      app.TimestepSpinner.Limits = [0 Inf];
      app.TimestepSpinner.Position = [104 140 100 22];
      app.TimestepSpinner.Value = 0.001;

      % Create DragEditFieldLabel
      app.DragEditFieldLabel = uilabel(app.LeftPanel);
      app.DragEditFieldLabel.HorizontalAlignment = 'right';
      app.DragEditFieldLabel.Position = [lmargin 110 31 22];
      app.DragEditFieldLabel.Text = 'Drag';

      % Create DragEditField
      app.DragEditField = uieditfield(app.LeftPanel, 'text');
      app.DragEditField.Position = [104 110 100 22];
      app.DragEditField.Value = '{eye(3), eye(3)}';

      % Create NumStepsSpinnerLabel
      app.NumStepsSpinnerLabel = uilabel(app.LeftPanel);
      app.NumStepsSpinnerLabel.HorizontalAlignment = 'right';
      app.NumStepsSpinnerLabel.Position = [lmargin 77 68 22];
      app.NumStepsSpinnerLabel.Text = 'Num. Steps';

      % Create NumStepsSpinner
      app.NumStepsSpinner = uispinner(app.LeftPanel);
      app.NumStepsSpinner.Position = [133 77 100 22];
      app.NumStepsSpinner.Value = 1000;

      % Create InitialPositionEditFieldLabel
      app.InitialPositionEditFieldLabel = uilabel(app.LeftPanel);
      app.InitialPositionEditFieldLabel.Position = [lmargin 323 79 22];
      app.InitialPositionEditFieldLabel.Text = 'Initial Position';

      % Create InitialPositionEditField
      app.InitialPositionEditField = uieditfield(app.LeftPanel, 'text');
      app.InitialPositionEditField.ValueChangedFcn = createCallbackFcn(app, @InitialPositionEditFieldValueChanged, true);
      app.InitialPositionEditField.Position = [132 323 100 22];
      app.InitialPositionEditField.Value = '[0, 0, 0]';

      % Create InitialRotationEditFieldLabel
      app.InitialRotationEditFieldLabel = uilabel(app.LeftPanel);
      app.InitialRotationEditFieldLabel.Position = [lmargin 291 82 22];
      app.InitialRotationEditFieldLabel.Text = 'Initial Rotation';

      % Create InitialRotationEditField
      app.InitialRotationEditField = uieditfield(app.LeftPanel, 'text');
      app.InitialRotationEditField.ValueChangedFcn = createCallbackFcn(app, @InitialRotationEditFieldValueChanged, true);
      app.InitialRotationEditField.Position = [132 291 100 22];
      app.InitialRotationEditField.Value = 'eye(3)';
    end
    
    function createRightComponents(app)
      
      lmargin = 10;
      width = 360;
      height = 180;

      % Create UIAxes
      app.UIAxes = uiaxes(app.RightPanel);
      title(app.UIAxes, 'Title')
      xlabel(app.UIAxes, 'X')
      ylabel(app.UIAxes, 'Y')
      app.UIAxes.Position = [lmargin height+15+30 width height];

      % Create UIAxes2
      app.UIAxes2 = uiaxes(app.RightPanel);
      xlabel(app.UIAxes2, 'X')
      ylabel(app.UIAxes2, 'Y')
      app.UIAxes2.Position = [lmargin 10 width height];

      % Create DropDown
      app.DropDown = uidropdown(app.RightPanel);
      app.DropDown.Items = {'Position', 'Rotation', 'Force', 'Torque'};
      app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
      app.DropDown.Position = [160 height+10 100 22];
      app.DropDown.Value = 'Position';
      
    end
  end
  
  methods (Access=public)
    function app=Isolated()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTwoColumn();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end