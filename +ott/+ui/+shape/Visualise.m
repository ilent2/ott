classdef Visualise < matlab.apps.AppBase
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
  end

  % Properties that correspond to app components
  properties (Access = public)
    UIFigure                    matlab.ui.Figure
    GridLayout                  matlab.ui.container.GridLayout
    LeftPanel                   matlab.ui.container.Panel
    ShapeDropDownLabel          matlab.ui.control.Label
    ShapeDropDown               matlab.ui.control.DropDown
    VisualisationDropDownLabel  matlab.ui.control.Label
    VisualisationDropDown       matlab.ui.control.DropDown
    UpdateButton                matlab.ui.control.Button
    AutoupdateCheckBox          matlab.ui.control.CheckBox
    RightPanel                  matlab.ui.container.Panel
    UIAxes                      matlab.ui.control.UIAxes
  end

  % Properties that correspond to apps with auto-reflow
  properties (Access = private)
    onePanelWidth = 576;
  end

  % Callbacks that handle component events
  methods (Access = private)

    % Changes arrangement of the app based on UIFigure width
    function updateAppLayout(app, event)
      
      warning('This should change to match CadFileLoader.m');
      
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {351, 351};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {220, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
    end
  end

  % Component initialization
  methods (Access = private)

    % Create UIFigure and components
    function createComponents(app)

      % Create UIFigure and hide until all components are created
      app.UIFigure = uifigure('Visible', 'off');
      app.UIFigure.AutoResizeChildren = 'off';
      app.UIFigure.Position = [100 100 634 351];
      app.UIFigure.Name = 'UI Figure';
      app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

      % Create GridLayout
      app.GridLayout = uigridlayout(app.UIFigure);
      app.GridLayout.ColumnWidth = {220, '1x'};
      app.GridLayout.RowHeight = {'1x'};
      app.GridLayout.ColumnSpacing = 0;
      app.GridLayout.RowSpacing = 0;
      app.GridLayout.Padding = [0 0 0 0];
      app.GridLayout.Scrollable = 'on';

      % Create LeftPanel
      app.LeftPanel = uipanel(app.GridLayout);
      app.LeftPanel.Layout.Row = 1;
      app.LeftPanel.Layout.Column = 1;

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

      % Create RightPanel
      app.RightPanel = uipanel(app.GridLayout);
      app.RightPanel.Layout.Row = 1;
      app.RightPanel.Layout.Column = 2;

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

      % Show the figure after all components are created
      app.UIFigure.Visible = 'on';
    end
  end

  % App creation and deletion
  methods (Access = public)

    function app = Visualise(shape)
      % Start the shape visualisation interface
      
      if nargin == 1
        % TODO
      end
      
      % Create UIFigure and components
      createComponents(app)

      % Register the app with App Designer
      registerApp(app, app.UIFigure)

      if nargout == 0
        clear app
      end
    end

    % Code that executes before app deletion
    function delete(app)

      % Delete UIFigure when app is deleted
      delete(app.UIFigure)
    end
  end
end