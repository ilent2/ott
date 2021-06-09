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

% TODO: Add support for input arguments

  properties (Constant)
    cnameText = 'ForcePosition';

    nameText = 'Generate Force-Position plot';

    aboutText = ['Generate plots of the force as a function of ' ...
      'displacement or rotation angle.'];
    
    helpText = {ott.ui.tools.ForcePosition.aboutText, ...
      '', ...
      'This interface assumes all units are SI base units, i.e. ', ...
      'meters for lengths, Newtons for forces.  Angles are in degrees.', ...
      '', ...
      'beam - Specify the beam variable name in the matlab workspace.', ...
      '', ...
      'tmatrix - Specify the T-matrix or Particle variable name in the', ...
      ' matlab workspace.', ...
      '', ...
      'direction - Direction to translate/rotate particle.', ...
      '', ...
      'resolution - Number of positions/rotations.', ...
      '', ...
      'range - Specifies distance to translate/rotate particle.  Uses', ...
      ' degrees for rotations and meters for translations.  Specify the', ...
      ' range as a start (left) and end (right) value.', ...
      '', ...
      'position - Initial position of particle.', ...
      '', ...
      'rotation - Initial rotation of particle.'};
    
    windowName = ott.ui.tools.ForcePosition.nameText;
    windowSize = [640, 400];
  end
  
  properties (Access=public)
    
    % Left panel
    LeftGrid                        matlab.ui.container.GridLayout
    BeamDropDown                    ott.ui.support.VariableDropDown
    TmatrixDropDown                 ott.ui.support.VariableDropDown
    DirectionDropDown               ott.ui.support.LabeledDropDown
    ResolutionSpinner               ott.ui.support.LabeledSpinner
    RangeSpinners                   ott.ui.support.RangeSpinners
    PositionXyzSpinner              ott.ui.support.XyzSpinners
    RotationXyzSpinner              ott.ui.support.XyzSpinners
    UpdateButton                    ott.ui.support.UpdateWithProgress
    
    % Right panel
    RightGrid                       matlab.ui.container.GridLayout
    ForceAxes                       matlab.ui.control.UIAxes
    TorqueAxes                      matlab.ui.control.UIAxes
    
  end
  
  methods (Access=protected)
    function progressCb(app, data)
      % Update the progress bar
      
      app.UpdateButton.Level = data.iteration/data.total * 100;
      drawnow;
    end
    
    function code = generateCode(app)
      code = {};
      code{end+1} = '% Required variables: beam, particle';
      code{end+1} = '';

      % Get range
      code{end+1} = ['theta = linspace(', num2str(app.RangeSpinners.Value(1)), ...
          ', ', num2str(app.RangeSpinners.Value(2)), ', ', ...
          num2str(app.ResolutionSpinner.Value), ');'];

      % Get initial position/rotation
      code{end+1} = ['xyz0 = [', num2str(app.PositionXyzSpinner.Value(:).'), '].'';'];
      code{end+1} = ['Rx = ', num2str(app.RotationXyzSpinner.Value(1)), ';'];
      code{end+1} = ['Ry = ', num2str(app.RotationXyzSpinner.Value(2)), ';'];
      code{end+1} = ['Rz = ', num2str(app.RotationXyzSpinner.Value(3)), ';'];
      code{end+1} = 'Rxyz0 = ott.utils.rotx(Rx)*ott.utils.roty(Ry)*ott.utils.rotz(Rz);';
      code{end+1} = '';

      switch app.DirectionDropDown.Value
        case 'x'
          code{end+1} = 'xyz = [1;0;0].*theta + xyz0;';
          code{end+1} = 'R = Rxyz0;';
          xtext = 'Position [m]';
        case 'y'
          code{end+1} = 'xyz = [0;1;0].*theta + xyz0;';
          code{end+1} = 'R = Rxyz0;';
          xtext = 'Position [m]';
        case 'z'
          code{end+1} = 'xyz = [0;0;1].*theta + xyz0;';
          code{end+1} = 'R = Rxyz0;';
          xtext = 'Position [m]';
        case 'Rx'
          code{end+1} = 'R = ott.utils.rotx(theta, ''usecell'', true);';
          code{end+1} = 'R = cellfun(@(r) r * R0, R, ''UniformOutput'', false);';
          code{end+1} = 'xyz = xyz0;';
          xtext = 'Angle [deg]';
        case 'Ry'
          code{end+1} = 'R = ott.utils.roty(theta, ''usecell'', true);';
          code{end+1} = 'R = cellfun(@(r) r * R0, R, ''UniformOutput'', false);';
          code{end+1} = 'xyz = xyz0;';
          xtext = 'Angle [deg]';
        case 'Rz'
          code{end+1} = 'R = ott.utils.rotz(theta, ''usecell'', true);';
          code{end+1} = 'R = cellfun(@(r) r * R0, R, ''UniformOutput'', false);';
          code{end+1} = 'xyz = xyz0;';
          xtext = 'Angle [deg]';
        otherwise
          error('Internal error');
      end
      code{end+1} = '';

      % Calculate force/torque
      code{end+1} = '[force, torque] = beam.forcetorque(particle, ...';
      code{end+1} = '    ''position'', xyz, ''rotation'', R);';
      code{end+1} = '';

      % Update plots
      code{end+1} = 'figure();';
      code{end+1} = 'subplot(1, 2, 1);';
      code{end+1} = 'plot(theta, -force);';
      code{end+1} = 'ylabel(''Force [N]'');';
      code{end+1} = ['xlabel(''', xtext, ''');'];
      code{end+1} = 'legend(''X'', ''Y'', ''Z'');';
      code{end+1} = 'subplot(1, 2, 2);';
      code{end+1} = 'plot(theta, -torque);';
      code{end+1} = 'ylabel(''Torque [Nm]'');';
      code{end+1} = ['xlabel(''', xtext, ''');'];
      code{end+1} = 'legend(''X'', ''Y'', ''Z'');';
    end
    
    function generatePlots(app)
      
      % Configure guage
      app.UpdateButton.Gauge.Enable = true;
      app.UpdateButton.Level = 0;
      app.UpdateButton.clearErrors();
      
      try
      
        % Get beam/particle
        tmatrix = app.TmatrixDropDown.Variable;
        beam = app.BeamDropDown.Variable;

        % Get range
        t0 = app.RangeSpinners.Value(1);
        t1 = app.RangeSpinners.Value(2);
        t = linspace(t0, t1, app.ResolutionSpinner.Value);

        % Get initial position/rotation
        xyz0 = app.PositionXyzSpinner.Value(:);
        Rxyz0 = app.RotationXyzSpinner.Value;
        R0 = ott.utils.rotx(Rxyz0(1))*ott.utils.roty(Rxyz0(2))*ott.utils.rotz(Rxyz0(3));

        switch app.DirectionDropDown.Value
          case 'x'
            xyz = [1;0;0].*t + xyz0;
            R = R0;
          case 'y'
            xyz = [0;1;0].*t + xyz0;
            R = R0;
          case 'z'
            xyz = [0;0;1].*t + xyz0;
            R = R0;
          case 'Rx'
            R = ott.utils.rotx(t, 'usecell', true);
            R = cellfun(@(r) r * R0, R, 'UniformOutput', false);
            xyz = xyz0;
          case 'Ry'
            R = ott.utils.roty(t, 'usecell', true);
            R = cellfun(@(r) r * R0, R, 'UniformOutput', false);
            xyz = xyz0;
          case 'Rz'
            R = ott.utils.rotz(t, 'usecell', true);
            R = cellfun(@(r) r * R0, R, 'UniformOutput', false);
            xyz = xyz0;
          otherwise
            error('Internal error');
        end

        % Calculate force/torque
        [force, torque] = beam.forcetorque(tmatrix, ...
            'position', xyz, 'rotation', R, ...
            'progress', @(a) app.progressCb(a));
      
        % Update plots
        plot(app.ForceAxes, t, -force);
        legend('X', 'Y', 'Z');
        plot(app.TorqueAxes, t, -torque);
        legend('X', 'Y', 'Z');
          
      catch ME
        app.UpdateButton.setError();
        rethrow(ME);
      end
      
      % Update guage
      app.UpdateButton.Gauge.Enable = false;
      
    end
    
    function setDefaultValues(app, ~)
      % Reset window to defaults
      
      % Set widget values
      app.BeamDropDown.Value = '';
      app.TmatrixDropDown.Value = '';
      app.DirectionDropDown.Value = 'z';
      app.ResolutionSpinner.Value = 100;
      % app.RangeSpinners.Value = [-1e-6, 1e-6]; % in directionValueChCb
      app.PositionXyzSpinner.Value = [0,0,0];
      app.RotationXyzSpinner.Value = [0,0,0];
      app.UpdateButton.Level = 0;
      app.UpdateButton.clearErrors();
      
      % Update direction specific things
      app.directionValueChangedCb();
      
      % Clear plots
      cla(app.ForceAxes);
      cla(app.TorqueAxes);
    end
    
    function directionValueChangedCb(app)

      % Update graph labels
      switch app.DirectionDropDown.Value
        case 'x'
          xlabel(app.ForceAxes, 'Position (m)')
          xlabel(app.TorqueAxes, 'Position (m)')
          app.RangeSpinners.Step = 5e-7;
          app.RangeSpinners.Value = [-1e-6, 1e-6];
        case 'y'
          xlabel(app.ForceAxes, 'Position (m)')
          xlabel(app.TorqueAxes, 'Position (m)')
          app.RangeSpinners.Step = 5e-7;
          app.RangeSpinners.Value = [-1e-6, 1e-6];
        case 'z'
          xlabel(app.ForceAxes, 'Position (m)')
          xlabel(app.TorqueAxes, 'Position (m)')
          app.RangeSpinners.Step = 5e-7;
          app.RangeSpinners.Value = [-1e-6, 1e-6];
        case 'Rx'
          xlabel(app.ForceAxes, 'Angle (deg)')
          xlabel(app.TorqueAxes, 'Angle (deg)')
          app.RangeSpinners.Step = 10;
          app.RangeSpinners.Value = [-90, 90];
        case 'Ry'
          xlabel(app.ForceAxes, 'Angle (deg)')
          xlabel(app.TorqueAxes, 'Angle (deg)')
          app.RangeSpinners.Step = 10;
          app.RangeSpinners.Value = [-90, 90];
        case 'Rz'
          xlabel(app.ForceAxes, 'Angle (deg)')
          xlabel(app.TorqueAxes, 'Angle (deg)')
          app.RangeSpinners.Step = 10;
          app.RangeSpinners.Value = [-90, 90];
        otherwise
          error('Internal error');
      end
    end
    
    function createLeftComponents(app)
      
      wwidth = 120;
      
      % Left layout grid
      app.LeftGrid = uigridlayout(app.LeftPanel);
      app.LeftGrid.RowHeight = repmat({32}, 1, 9);
      app.LeftGrid.RowHeight{end-1} = '1x';
      app.LeftGrid.ColumnWidth = {'1x'};
      app.LeftGrid.RowSpacing = 0;
      
      % Beam selection
      app.BeamDropDown = ott.ui.support.VariableDropDown(app.LeftGrid, ...
        'label', 'Beam', 'wwidth', wwidth, 'filter', 'ott.beam.Beam');
      app.BeamDropDown.Layout.Row = 1;
      app.BeamDropDown.Layout.Column = 1;
      app.registerRefreshInput(app.BeamDropDown);
      
      % T-matrix selection
      app.TmatrixDropDown = ott.ui.support.VariableDropDown(app.LeftGrid, ...
        'label', 'Particle', 'wwidth', wwidth, ...
        'filter', {'ott.tmatrix.Tmatrix', 'ott.particle.Particle'});
      app.TmatrixDropDown.Layout.Row = 2;
      app.TmatrixDropDown.Layout.Column = 1;
      app.registerRefreshInput(app.TmatrixDropDown);
      
      % Direction
      app.DirectionDropDown = ott.ui.support.LabeledDropDown(app.LeftGrid, ...
        'label', 'Direction', 'wwidth', wwidth);
      app.DirectionDropDown.Items = {'X Translation', 'Y Translation', ...
        'Z Translation', 'X Rotation', 'Y Rotation', 'Z Rotation'};
      app.DirectionDropDown.ItemsData = {'x', 'y', 'z', 'Rx', 'Ry', 'Rz'};
      app.DirectionDropDown.Layout.Row = 3;
      app.DirectionDropDown.Layout.Column = 1;
      app.DirectionDropDown.ValueChangedFcn = @(~,~) app.directionValueChangedCb();
      
      % Resolution
      app.ResolutionSpinner = ott.ui.support.LabeledSpinner(app.LeftGrid, ...
        'label', 'Resolution', 'wwidth', wwidth);
      app.ResolutionSpinner.Step = 1;
      app.ResolutionSpinner.Limits = [1, Inf];
      app.ResolutionSpinner.Layout.Row = 4;
      app.ResolutionSpinner.Layout.Column = 1;
      app.ResolutionSpinner.RoundFractionalValues = true;
      
      % Range
      app.RangeSpinners = ott.ui.support.RangeSpinners(app.LeftGrid, ...
          'wwidth', 160);
      app.RangeSpinners.Step = 5e-7;
      app.RangeSpinners.Layout.Row = 5;
      app.RangeSpinners.Layout.Column = 1;

      % Create Position
      app.PositionXyzSpinner = ott.ui.support.XyzSpinners(app.LeftGrid, ...
        'label', 'Position');
      app.PositionXyzSpinner.Layout.Row = 6;
      app.PositionXyzSpinner.Layout.Column = 1;
      app.PositionXyzSpinner.Step = 5e-7;

      % Create Rotation
      app.RotationXyzSpinner = ott.ui.support.XyzSpinners(app.LeftGrid, ...
          'label', 'Rotation');
      app.RotationXyzSpinner.Layout.Row = 7;
      app.RotationXyzSpinner.Layout.Column = 1;
      app.RotationXyzSpinner.Step = 10;

      % Create CalculateButton
      app.UpdateButton = ott.ui.support.UpdateWithProgress(app.LeftGrid);
      app.UpdateButton.Layout.Row = 9;
      app.UpdateButton.Layout.Column = 1;
      app.UpdateButton.UpdateCalledFcn = @(~,~) app.generatePlots();
    end
    
    function createRightComponents(app)
      
      % Left layout grid
      app.RightGrid = uigridlayout(app.RightPanel);
      app.RightGrid.RowHeight = {'1x', '1x'};
      app.RightGrid.ColumnWidth = {'1x'};
      app.RightGrid.RowSpacing = 20;
      
      % Create UIAxes
      app.ForceAxes = uiaxes(app.RightGrid);
      xlabel(app.ForceAxes, 'Position (m)')
      ylabel(app.ForceAxes, 'Force (N)')
      app.ForceAxes.Layout.Row = 1;
      app.ForceAxes.Layout.Column = 1;

      % Create UIAxes2
      app.TorqueAxes = uiaxes(app.RightGrid);
      xlabel(app.TorqueAxes, 'Position (m)')
      ylabel(app.TorqueAxes, 'Torque (Nm)')
      app.TorqueAxes.Layout.Row = 2;
      app.TorqueAxes.Layout.Column = 1;
      
    end
  end
  
  methods (Access=public)
    function app = ForcePosition()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTwoColumn();
      
      if nargout == 0
        clear app;
      end
    end
  end
end
