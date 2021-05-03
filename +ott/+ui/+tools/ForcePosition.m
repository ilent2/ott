classdef ForcePosition < ott.ui.support.AppTwoColumn ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Generate a plot of the force as a function of particle position.
%
% This app can be launched from the Launcher via Tools -> ForcePosition
% or directly from the command line with:
%
%   ott.ui.tools.ForcePosition()

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'ForcePosition';

    nameText = 'Generate Force-Position plot';

    aboutText = ['Generate plots of the force as a function of ' ...
      'displacement or rotation angle.'];
    
    helpText = {ott.ui.tools.ForcePosition.aboutText, ...
      ''};
    
    windowName = ott.ui.tools.ForcePosition.nameText;
    windowSize = [640, 460];
  end
  
  properties (Access=public)
    
    % Left panel
    LeftGrid                        matlab.ui.control.GridLayout
    BeamDropDown                    ott.ui.support.VariableDropdown
    TmatrixDropDown                 ott.ui.support.VariableDropdown
    DirectionDropdown               ott.ui.support.LabeledDropdown
    ResolutionSpinner               ott.ui.support.LabeledSpinner
    RangeSpinners                   ott.ui.support.RangeSpinners
    PositionXyzSpinner              ott.ui.support.XyzSpinners
    RotationXyzSpinner              ott.ui.support.XyzSpinners
    UpdateButton                    ott.ui.support.UpdateWithProgress
    
    % Right panel
    RightGrid                       matlab.ui.control.GridLayout
    ForceAxes                       matlab.ui.control.UIAxes
    TorqueAxes                      matlab.ui.control.UIAxes
    
  end
  
  methods (Access=protected)
    function createLeftComponents(app)

      % Create CalculateButton
      app.CalculateButton = uibutton(app.LeftPanel, 'push');
      app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
      app.CalculateButton.Position = [148 16 71 22];
      app.CalculateButton.Text = 'Calculate';

      % Create DirectionDropDownLabel
      app.DirectionDropDownLabel = uilabel(app.LeftPanel);
      app.DirectionDropDownLabel.HorizontalAlignment = 'right';
      app.DirectionDropDownLabel.Position = [9 296 53 22];
      app.DirectionDropDownLabel.Text = 'Direction';

      % Create DirectionDropDown
      app.DirectionDropDown = uidropdown(app.LeftPanel);
      app.DirectionDropDown.Items = {'X Translation', 'Y Translation', 'Z Translation', 'X Rotation', 'Y Rotation', 'Z Rotation'};
      app.DirectionDropDown.ItemsData = {'x', 'y', 'z', 'Rx', 'Ry', 'Rz'};
      app.DirectionDropDown.Position = [108 296 100 22];
      app.DirectionDropDown.Value = 'z';

      % Create ResolutionSpinnerLabel
      app.ResolutionSpinnerLabel = uilabel(app.LeftPanel);
      app.ResolutionSpinnerLabel.HorizontalAlignment = 'right';
      app.ResolutionSpinnerLabel.Position = [9 222 62 22];
      app.ResolutionSpinnerLabel.Text = 'Resolution';

      % Create ResolutionSpinner
      app.ResolutionSpinner = uispinner(app.LeftPanel);
      app.ResolutionSpinner.Position = [108 222 100 22];
      app.ResolutionSpinner.Value = 100;

      % Create RangeSpinnerLabel
      app.RangeSpinnerLabel = uilabel(app.LeftPanel);
      app.RangeSpinnerLabel.HorizontalAlignment = 'right';
      app.RangeSpinnerLabel.Position = [9 259 41 22];
      app.RangeSpinnerLabel.Text = 'Range';

      % Create RangeSpinner
      app.RangeSpinner = uispinner(app.LeftPanel);
      app.RangeSpinner.Position = [65 259 67 22];
      app.RangeSpinner.Value = -2;

      % Create RangeSpinner_2
      app.RangeSpinner_2 = uispinner(app.LeftPanel);
      app.RangeSpinner_2.Position = [139 259 69 22];
      app.RangeSpinner_2.Value = 2;

      % Create Gauge
      app.Gauge = uigauge(app.LeftPanel, 'linear');
      app.Gauge.Enable = 'off';
      app.Gauge.FontSize = 8;
      app.Gauge.Position = [9 13 133 29];

      % Create BeamDropDownLabel
      app.BeamDropDownLabel = uilabel(app.LeftPanel);
      app.BeamDropDownLabel.HorizontalAlignment = 'right';
      app.BeamDropDownLabel.Position = [9 417 37 22];
      app.BeamDropDownLabel.Text = 'Beam';

      % Create BeamDropDown
      app.BeamDropDown = uidropdown(app.LeftPanel);
      app.BeamDropDown.Items = {};
      app.BeamDropDown.Editable = 'on';
      app.BeamDropDown.BackgroundColor = [1 1 1];
      app.BeamDropDown.Position = [108 417 100 22];
      app.BeamDropDown.Value = {};

      % Create TmatrixDropDownLabel
      app.TmatrixDropDownLabel = uilabel(app.LeftPanel);
      app.TmatrixDropDownLabel.HorizontalAlignment = 'right';
      app.TmatrixDropDownLabel.Position = [9 380 49 22];
      app.TmatrixDropDownLabel.Text = 'T-matrix';

      % Create TmatrixDropDown
      app.TmatrixDropDown = uidropdown(app.LeftPanel);
      app.TmatrixDropDown.Items = {};
      app.TmatrixDropDown.Editable = 'on';
      app.TmatrixDropDown.BackgroundColor = [1 1 1];
      app.TmatrixDropDown.Position = [108 380 100 22];
      app.TmatrixDropDown.Value = {};

      % Create InitialPositionEditFieldLabel
      app.InitialPositionEditFieldLabel = uilabel(app.LeftPanel);
      app.InitialPositionEditFieldLabel.HorizontalAlignment = 'right';
      app.InitialPositionEditFieldLabel.Position = [9 175 80 22];
      app.InitialPositionEditFieldLabel.Text = 'Initial Position';

      % Create InitialPositionEditField
      app.InitialPositionEditField = uieditfield(app.LeftPanel, 'text');
      app.InitialPositionEditField.Position = [108 175 100 22];
      app.InitialPositionEditField.Value = '[0;0;0]';

      % Create InitialRotationEditFieldLabel
      app.InitialRotationEditFieldLabel = uilabel(app.LeftPanel);
      app.InitialRotationEditFieldLabel.HorizontalAlignment = 'right';
      app.InitialRotationEditFieldLabel.Position = [9 138 82 22];
      app.InitialRotationEditFieldLabel.Text = 'Initial Rotation';

      % Create InitialRotationEditField
      app.InitialRotationEditField = uieditfield(app.LeftPanel, 'text');
      app.InitialRotationEditField.Position = [108 138 100 22];
      app.InitialRotationEditField.Value = 'eye(3)';
    end
    
    function createRightComponents(app)
      
      width = 360;
      height = 220;
      
      % Create UIAxes
      app.ForceAxes = uiaxes(app.RightPanel);
      xlabel(app.ForceAxes, 'X')
      ylabel(app.ForceAxes, 'Force')
      app.ForceAxes.Position = [10 height+10 width height];

      % Create UIAxes2
      app.TorqueAxes = uiaxes(app.RightPanel);
      xlabel(app.TorqueAxes, 'X')
      ylabel(app.TorqueAxes, 'Torque')
      app.TorqueAxes.Position = [10 10 width height];
    end
  end
  
  methods (Access=public)
    function app=ForcePosition()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTwoColumn();
      
      if nargout == 0
        clear app;
      end
    end
  end
end
