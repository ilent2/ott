classdef Isolated < ott.ui.support.AppTwoColumn ...
    & ott.ui.support.RefreshInputsMenu ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.AppProducer
% Simulate the dynamics of an isolated particle.
% Internally uses the ott.dynamics.Isolated class.
%
% This application can be run from the Launcher via Dynamics -> Isolated
% or launched directly from the command line with:
%
%   ott.ui.dynamics.Isolated()

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Isolated';

    nameText = 'Simulate Isolated Dynamics';

    aboutText = ['Simulate dynamics of an isolated particle.'];
    
    helpText = {ott.ui.dynamics.Isolated.aboutText, ...
      ''};
    
    windowName = ott.ui.dynamics.Isolated.nameText;
    windowSize = [640, 420];
  end
  
  properties (Access=public)
    % Left panel
    LeftGrid            matlab.ui.container.GridLayout
    BeamDropDown        ott.ui.support.VariableDropDown
    ParticleDropDown    ott.ui.support.VariableDropDown
    TimestepSpinner     ott.ui.support.LabeledSpinner
    NumStepsSpinner     ott.ui.support.LabeledSpinner
    PositionXyzSpinner  ott.ui.support.XyzSpinners
    RotationXyzSpinner  ott.ui.support.XyzSpinners
    
    % Right panel
    RightGrid         matlab.ui.container.GridLayout
    TopAxes           matlab.ui.control.UIAxes
    BottomAxes        matlab.ui.control.UIAxes
    PlotType          matlab.ui.control.DropDown
  end
  
  methods (Access=protected)
    function setDefaultValues(app)
      
      % Left panel
      app.VariableName.Value = '';
      app.BeamDropDown.Value = '';
      app.ParticleDropDown.Value = '';
      app.TimestepSpinner.Value = 1e-3;
      app.NumStepsSpinner.Value = 100;
      app.PositionXyzSpinner.Value = [0,0,0];
      app.RotationXyzSpinner.Value = [0,0,0];
      app.UpdateButton.Level = 0;
      
      % Right panel
      app.PlotType.Value = 'Position';
    end
    
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function data = generateData(app)
      
      % Set progress bar to 0
      
      % Start simulation
      % Update progress bar every N seconds
      
      data = []; % TODO
      
    end
    
    function createLeftComponents(app)
      
      % Create grid
      app.LeftGrid = uigridlayout(app.LeftPanel);
      app.LeftGrid.ColumnWidth = {'1x'};
      app.LeftGrid.RowHeight = repmat({32}, 1, 9);
      app.LeftGrid.RowHeight{end-1} = '1x';
      
      % Variable name
      app.VariableName = ott.ui.support.OutputVariableEntry(app.LeftGrid);
      app.VariableName.Layout.Row = 1;
      app.VariableName.Layout.Column = 1;
      
      % Beam selector
      app.BeamDropDown = ott.ui.support.VariableDropDown(...
        app.LeftGrid, 'label', 'Beam');
      app.BeamDropDown.Layout.Row = 2;
      app.BeamDropDown.Layout.Column = 1;
      app.BeamDropDown.ValueChangedFcn = createCallbackFcn(app, ...
          @UpdateParametersCb, true);
      
      % Particle selector
      app.ParticleDropDown = ott.ui.support.VariableDropDown(...
        app.LeftGrid, 'label', 'Particle');
      app.ParticleDropDown.Layout.Row = 3;
      app.ParticleDropDown.Layout.Column = 1;
      app.ParticleDropDown.ValueChangedFcn = createCallbackFcn(app, ...
          @UpdateParametersCb, true);
      
      % Time step size
      app.TimestepSpinner = ott.ui.support.LabeledSpinner(app.LeftGrid, ...
          'label', 'Timestep');
      app.TimestepSpinner.Layout.Row = 4;
      app.TimestepSpinner.Layout.Column = 1;
      app.TimestepSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @UpdateParametersCb, true);
      
      % Num time steps
      app.NumStepsSpinner = ott.ui.support.LabeledSpinner(app.LeftGrid, ...
          'label', 'Num. Steps');
      app.NumStepsSpinner.Layout.Row = 5;
      app.NumStepsSpinner.Layout.Column = 1;
      app.NumStepsSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @UpdateParametersCb, true);
      
      % Initial position
      app.PositionXyzSpinner = ott.ui.support.XyzSpinners(...
          app.LeftGrid, 'label', 'Position');
      app.PositionXyzSpinner.Layout.Row = 6;
      app.PositionXyzSpinner.Layout.Column = 1;
      app.PositionXyzSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @UpdateParametersCb, true);
      
      % Initial rotation
      app.RotationXyzSpinner = ott.ui.support.XyzSpinners(...
          app.LeftGrid, 'label', 'Rotation');
      app.RotationXyzSpinner.Layout.Row = 7;
      app.RotationXyzSpinner.Layout.Column = 1;
      app.RotationXyzSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @UpdateParametersCb, true);
    
      % Progress bar and update button
      app.UpdateButton = ott.ui.support.UpdateWithProgress(app.LeftGrid);
      app.UpdateButton.Layout.Row = 9;
      app.UpdateButton.Layout.Column = 1;
    end
    
    function createRightComponents(app)
      
      % Create grid
      app.RightGrid = uigridlayout(app.RightPanel);
      app.RightGrid.ColumnWidth = {'1x'};
      app.RightGrid.RowHeight = {'1x', 32, '1x'};

      % Create UIAxes
      app.TopAxes = uiaxes(app.RightGrid);
      title(app.TopAxes, 'Title');
      xlabel(app.TopAxes, 'X');
      ylabel(app.TopAxes, 'Y');
      app.TopAxes.Layout.Row = 1;
      app.TopAxes.Layout.Column = 1;
      
      rowGrid = uigridlayout(app.RightGrid);
      rowGrid.Padding = [0,0,0,10];
      rowGrid.ColumnWidth = {'1.2x', 100, '1x'};
      rowGrid.RowHeight = {22};
      rowGrid.Layout.Row = 2;
      rowGrid.Layout.Column = 1;

      % Create DropDown
      app.PlotType = uidropdown(rowGrid);
      app.PlotType.Items = {'Position', 'Rotation', 'Force', 'Torque'};
      app.PlotType.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
      app.PlotType.Layout.Row = 1;
      app.PlotType.Layout.Column = 2;

      % Create UIAxes2
      app.BottomAxes = uiaxes(app.RightGrid);
      xlabel(app.BottomAxes, 'X');
      ylabel(app.BottomAxes, 'Y');
      app.BottomAxes.Layout.Row = 3;
      app.BottomAxes.Layout.Column = 1;
      
    end
  end
  
  methods (Access=public)
    function app=Isolated()
      % Start the ForcePosition GUI
      
      % Construct widgets first
      app = app@ott.ui.support.AppTwoColumn();
      
      % Then connect producer callbacks
      app = app@ott.ui.support.AppProducer();
      
      if nargout == 0
        clear app;
      end
    end
  end
end