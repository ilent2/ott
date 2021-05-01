classdef Scattered < ott.ui.beam.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Scattered';

    nameText = 'Scatter Beam';

    aboutText = ['Generate a scattered beam from a particle and incident' ...
      'beam instance.'];
    
    helpText = {ott.ui.beam.Scattered.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 360];
  end
  
  properties (Access=public)
    GenerateButton               matlab.ui.control.Button
    IncidentBeamDropDownLabel    matlab.ui.control.Label
    IncidentBeamDropDown         matlab.ui.control.DropDown
    TmatrixDropDownLabel         matlab.ui.control.Label
    TmatrixDropDown              matlab.ui.control.DropDown
    ScatteredBeamEditFieldLabel  matlab.ui.control.Label
    ScatteredBeamEditField       matlab.ui.control.EditField
    PositionEditFieldLabel       matlab.ui.control.Label
    PositionEditField            matlab.ui.control.EditField
    RotationEditFieldLabel       matlab.ui.control.Label
    RotationEditField            matlab.ui.control.EditField
  end
  
  methods (Access=protected)
    function startupFcn(app)
    end
    
    function createRightComponents(app)
      app.createBeamPreview();
    end
    
    function createLeftComponents(app)

      % Create GenerateButton
      app.GenerateButton = uibutton(app.LeftPanel, 'push');
      app.GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateButtonPushed, true);
      app.GenerateButton.Position = [116 12 100 22];
      app.GenerateButton.Text = 'Generate';

      % Create IncidentBeamDropDownLabel
      app.IncidentBeamDropDownLabel = uilabel(app.LeftPanel);
      app.IncidentBeamDropDownLabel.HorizontalAlignment = 'right';
      app.IncidentBeamDropDownLabel.Position = [11 227 82 22];
      app.IncidentBeamDropDownLabel.Text = 'Incident Beam';

      % Create IncidentBeamDropDown
      app.IncidentBeamDropDown = uidropdown(app.LeftPanel);
      app.IncidentBeamDropDown.Items = {};
      app.IncidentBeamDropDown.Editable = 'on';
      app.IncidentBeamDropDown.BackgroundColor = [1 1 1];
      app.IncidentBeamDropDown.Position = [116 227 100 22];
      app.IncidentBeamDropDown.Value = {};

      % Create TmatrixDropDownLabel
      app.TmatrixDropDownLabel = uilabel(app.LeftPanel);
      app.TmatrixDropDownLabel.HorizontalAlignment = 'right';
      app.TmatrixDropDownLabel.Position = [11 190 49 22];
      app.TmatrixDropDownLabel.Text = 'T-matrix';

      % Create TmatrixDropDown
      app.TmatrixDropDown = uidropdown(app.LeftPanel);
      app.TmatrixDropDown.Items = {};
      app.TmatrixDropDown.Editable = 'on';
      app.TmatrixDropDown.BackgroundColor = [1 1 1];
      app.TmatrixDropDown.Position = [116 190 100 22];
      app.TmatrixDropDown.Value = {};

      % Create ScatteredBeamEditFieldLabel
      app.ScatteredBeamEditFieldLabel = uilabel(app.LeftPanel);
      app.ScatteredBeamEditFieldLabel.HorizontalAlignment = 'right';
      app.ScatteredBeamEditFieldLabel.Position = [11 143 92 22];
      app.ScatteredBeamEditFieldLabel.Text = 'Scattered Beam';

      % Create ScatteredBeamEditField
      app.ScatteredBeamEditField = uieditfield(app.LeftPanel, 'text');
      app.ScatteredBeamEditField.Position = [116 143 100 22];
      app.ScatteredBeamEditField.Value = 'ScatteredBeam';

      % Create PositionEditFieldLabel
      app.PositionEditFieldLabel = uilabel(app.LeftPanel);
      app.PositionEditFieldLabel.HorizontalAlignment = 'right';
      app.PositionEditFieldLabel.Position = [11 96 48 22];
      app.PositionEditFieldLabel.Text = 'Position';

      % Create PositionEditField
      app.PositionEditField = uieditfield(app.LeftPanel, 'text');
      app.PositionEditField.Position = [116 96 100 22];
      app.PositionEditField.Value = '[0, 0, 0]';

      % Create RotationEditFieldLabel
      app.RotationEditFieldLabel = uilabel(app.LeftPanel);
      app.RotationEditFieldLabel.HorizontalAlignment = 'right';
      app.RotationEditFieldLabel.Position = [11 59 50 22];
      app.RotationEditFieldLabel.Text = 'Rotation';

      % Create RotationEditField
      app.RotationEditField = uieditfield(app.LeftPanel, 'text');
      app.RotationEditField.Position = [116 59 100 22];
      app.RotationEditField.Value = 'eye(3)';
    end
  end
  
  methods (Access=public)
    function app=Scattered()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.AppBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end