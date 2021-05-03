classdef Visualise < ott.ui.support.AppTwoColumn
% Generate a visualisation of a OTT geometric shape.
%
% This GUI can be launched from the launcher under
% Shape -> Visualise or running the following command:
%
%   ott.ui.shape.Visualise()
%
% Or with a specific shape:
%
%   ott.ui.shape.Visualise(shape)

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Visualise';

    nameText = 'Visualise Shape';

    aboutText = ['Generate a visualisation of a shape.'];
    
    helpText = {ott.ui.shape.Visualise.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 350];
  end

  % Properties that correspond to app components
  properties (Access = public)
    ShapeDropDownLabel          matlab.ui.control.Label
    ShapeDropDown               matlab.ui.control.DropDown
    VisualisationDropDownLabel  matlab.ui.control.Label
    VisualisationDropDown       matlab.ui.control.DropDown
    UpdateButton                matlab.ui.control.Button
    AutoupdateCheckBox          matlab.ui.control.CheckBox
    UIAxes                      matlab.ui.control.UIAxes
  end

  % Properties that correspond to apps with auto-reflow
  properties (Access = private)
    onePanelWidth = 576;
  end

  % Component initialization
  methods (Access = protected)
    
    function startupFcn(app)
    end

    % Create UIFigure and components
    function createLeftComponents(app)

      % Create ShapeDropDownLabel
      app.ShapeDropDownLabel = uilabel(app.LeftPanel);
      app.ShapeDropDownLabel.HorizontalAlignment = 'right';
      app.ShapeDropDownLabel.Position = [45 312 40 22];
      app.ShapeDropDownLabel.Text = 'Shape';

      % Create ShapeDropDown
      app.ShapeDropDown = uidropdown(app.LeftPanel);
      app.ShapeDropDown.Editable = 'on';
      app.ShapeDropDown.BackgroundColor = [1 1 1];
      app.ShapeDropDown.Position = [100 312 100 22];

      % Create VisualisationDropDownLabel
      app.VisualisationDropDownLabel = uilabel(app.LeftPanel);
      app.VisualisationDropDownLabel.HorizontalAlignment = 'right';
      app.VisualisationDropDownLabel.Position = [13 277 72 22];
      app.VisualisationDropDownLabel.Text = 'Visualisation';

      % Create VisualisationDropDown
      app.VisualisationDropDown = uidropdown(app.LeftPanel);
      app.VisualisationDropDown.Items = {'Voxels', 'Surface'};
      app.VisualisationDropDown.Position = [100 277 100 22];
      app.VisualisationDropDown.Value = 'Surface';

      % Create UpdateButton
      app.UpdateButton = uibutton(app.LeftPanel, 'push');
      app.UpdateButton.Enable = 'off';
      app.UpdateButton.Position = [116 212 84 22];
      app.UpdateButton.Text = 'Update';

      % Create AutoupdateCheckBox
      app.AutoupdateCheckBox = uicheckbox(app.LeftPanel);
      app.AutoupdateCheckBox.Text = 'Auto-update';
      app.AutoupdateCheckBox.Position = [13 212 87 22];
      app.AutoupdateCheckBox.Value = true;
      
    end
    
    function createRightComponents(app)

      % Create UIAxes
      app.UIAxes = uiaxes(app.RightPanel);
      title(app.UIAxes, 'Preview')
      xlabel(app.UIAxes, '')
      ylabel(app.UIAxes, '')
      app.UIAxes.XAxisLocation = 'origin';
      app.UIAxes.XTick = [];
      app.UIAxes.YAxisLocation = 'origin';
      app.UIAxes.YTick = [];
      app.UIAxes.Position = [20 6 373 328];
    end
  end

  % App creation and deletion
  methods (Access = public)

    function app = Visualise(shape)
      % Start the shape visualisation interface
      
      app = app@ott.ui.shape.AppBase();
      
      if nargin == 1
        % TODO
      end
      
      if nargout == 0
        clear app;
      end
    end
  end
end