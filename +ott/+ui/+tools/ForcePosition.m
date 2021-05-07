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

% TODO: Change range step when direction changes
% TODO: Add support for input arguments

  properties (Constant)
    cnameText = 'ForcePosition';

    nameText = 'Generate Force-Position plot';

    aboutText = ['Generate plots of the force as a function of ' ...
      'displacement or rotation angle.'];
    
    helpText = {ott.ui.tools.ForcePosition.aboutText, ...
      ''};
    
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
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function setDefaultValues(app, ~)
      app.BeamDropDown.Value = '';
      app.TmatrixDropDown.Value = '';
      app.DirectionDropDown.Value = 'z';
      app.ResolutionSpinner.Value = 100;
      app.RangeSpinners.Value = [-1e-6, 1e-6];
      app.PositionXyzSpinner.Value = [0,0,0];
      app.RotationXyzSpinner.Value = [0,0,0];
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
        'label', 'Beam', 'wwidth', wwidth);
      app.BeamDropDown.Layout.Row = 1;
      app.BeamDropDown.Layout.Column = 1;
      app.registerRefreshInput(app.BeamDropDown);
      
      % T-matrix selection
      app.TmatrixDropDown = ott.ui.support.VariableDropDown(app.LeftGrid, ...
        'label', 'Particle', 'wwidth', wwidth);
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
      
      % Resolution
      app.ResolutionSpinner = ott.ui.support.LabeledSpinner(app.LeftGrid, ...
        'label', 'Resolution', 'wwidth', wwidth);
      app.ResolutionSpinner.Step = 1;
      app.ResolutionSpinner.Limits = [1, Inf];
      app.ResolutionSpinner.Layout.Row = 4;
      app.ResolutionSpinner.Layout.Column = 1;
      
      % Range
      app.RangeSpinners = ott.ui.support.RangeSpinners(app.LeftGrid);
      app.RangeSpinners.Step = 1e-7;
      app.RangeSpinners.Layout.Row = 5;
      app.RangeSpinners.Layout.Column = 1;

      % Create Position
      app.PositionXyzSpinner = ott.ui.support.XyzSpinners(app.LeftGrid, ...
        'label', 'Position');
      app.PositionXyzSpinner.Layout.Row = 6;
      app.PositionXyzSpinner.Layout.Column = 1;

      % Create Rotation
      app.RotationXyzSpinner = ott.ui.support.XyzSpinners(app.LeftGrid, ...
          'label', 'Rotation');
      app.RotationXyzSpinner.Layout.Row = 7;
      app.RotationXyzSpinner.Layout.Column = 1;

      % Create CalculateButton
      app.UpdateButton = ott.ui.support.UpdateWithProgress(app.LeftGrid);
      app.UpdateButton.Layout.Row = 9;
      app.UpdateButton.Layout.Column = 1;
    end
    
    function createRightComponents(app)
      
      % Left layout grid
      app.RightGrid = uigridlayout(app.RightPanel);
      app.RightGrid.RowHeight = {'1x', '1x'};
      app.RightGrid.ColumnWidth = {'1x'};
      app.RightGrid.RowSpacing = 20;
      
      % Create UIAxes
      app.ForceAxes = uiaxes(app.RightGrid);
      xlabel(app.ForceAxes, 'X')
      ylabel(app.ForceAxes, 'Force')
      app.ForceAxes.Layout.Row = 1;
      app.ForceAxes.Layout.Column = 1;

      % Create UIAxes2
      app.TorqueAxes = uiaxes(app.RightGrid);
      xlabel(app.TorqueAxes, 'X')
      ylabel(app.TorqueAxes, 'Torque')
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
