classdef CadFileLoader < ott.ui.shape.AppBase
% Generate a OTT shape from a computer aided design (CAD) file.
%
% This GUI can be launched from the launcher under
% Shape -> CadFileLoader or running the following command:
%
%   ott.ui.shape.CadFileLoader()
%
% The shape is stored internally and/or written to the matlab workspace
% if a variable name is given for the shape.  To access the internal
% instance use:
%
%   app = ott.ui.shape.CadFileLoader()
%   shape = app.shape

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'CadFileLoader';

    nameText = 'CAD File Loader';

    aboutText = ['Create an OTT shape from a computer aided design (CAD)', ...
      ' file.  Supports loading STL and Wavefront OBJ files.  Free tools', ...
      ' such as Blender can be used to create and convert between CAD', ...
      ' files.'];
    
    helpText = {ott.ui.shape.CadFileLoader.aboutText, ...
      '', ...
      ['VariableName (optional) : Name for the variable to ' ...
      'genrate in the Matlab workspace.'], ...
      '', ...
      'CAD File : Input CAD file.  Must be STL or OBJ file.', ...
      '', ...
      'Scale : Object scale factor.', ...
      '', ...
      'Offset (x, y, z) : Object offset, meters.', ...
      '', ...
      'Rotation (x, y, z) : Rotation, degrees.', ...
      '', ...
      ['Rotational/XY Symmetry: Set these if the objects has', ...
      'symetries (including rotations/translations).  If in ', ...
      'doubt, do not set these.']};
  end
  
  % Properties that correspond to app components
  properties (Access = public)
    UIFigure                        matlab.ui.Figure
    GridLayout                      matlab.ui.container.GridLayout
    LeftPanel                       matlab.ui.container.Panel
    CADFileEditFieldLabel           matlab.ui.control.Label
    CADFileEditField                matlab.ui.control.EditField
    Button                          matlab.ui.control.Button
    ScaleSpinnerLabel               matlab.ui.control.Label
    ScaleSpinner                    matlab.ui.control.Spinner
    XYMirrorSymmetryCheckBox        matlab.ui.control.CheckBox
    RotationalSymmetrySpinnerLabel  matlab.ui.control.Label
    RotationalSymmetrySpinner       matlab.ui.control.Spinner
    VariableNameEditFieldLabel      matlab.ui.control.Label
    VariableNameEditField           matlab.ui.control.EditField
    ShowPreviewCheckBox             matlab.ui.control.CheckBox
    UpdateButton                    matlab.ui.control.Button
    OffsetSpinnerLabel              matlab.ui.control.Label
    OffsetSpinner                   matlab.ui.control.Spinner
    OffsetSpinner_2                 matlab.ui.control.Spinner
    OffsetSpinner_3                 matlab.ui.control.Spinner
    OffsetSpinner_4                 matlab.ui.control.Spinner
    OffsetSpinner_5                 matlab.ui.control.Spinner
    RotationSpinnerLabel            matlab.ui.control.Label
    RotationSpinner                 matlab.ui.control.Spinner
    AutoupdateCheckBox              matlab.ui.control.CheckBox
    RightPanel                      matlab.ui.container.Panel
    UIAxes                          matlab.ui.control.UIAxes
  end
  
  % Properties that correspond to apps with auto-reflow
  properties (Access = private)
    onePanelWidth = 576;
    onePanelHeight = 300;
  end

  % Callbacks that handle component events
  methods (Access = private)

    % Changes arrangement of the app based on UIFigure width
    function updateAppLayout(app, ~)
      currentFigureWidth = app.UIFigure.Position(3);
      currentFigureHeight = app.UIFigure.Position(4);
      if(currentFigureWidth <= app.onePanelWidth ...
          || currentFigureHeight <= app.onePanelHeight)
        % Change to a 2x1 grid
        app.GridLayout.RowHeight = {418, 418};
        app.GridLayout.ColumnWidth = {'1x'};
        app.RightPanel.Layout.Row = 2;
        app.RightPanel.Layout.Column = 1;
      else
        % Change to a 1x2 grid
        app.GridLayout.RowHeight = {'1x'};
        app.GridLayout.ColumnWidth = {251, '1x'};
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
      app.UIFigure.Position = [ott.ui.support.defaultXy, 640, 420];
      app.UIFigure.Name = app.nameText;
      app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);
      
      % Create menus
      ott.ui.support.createLauncherMenuItem(app.UIFigure);
      ott.ui.support.createHelpMenu(app.UIFigure, app.helpText);

      % Create GridLayout
      app.GridLayout = uigridlayout(app.UIFigure);
      app.GridLayout.ColumnWidth = {251, '1x'};
      app.GridLayout.RowHeight = {'1x'};
      app.GridLayout.ColumnSpacing = 0;
      app.GridLayout.RowSpacing = 0;
      app.GridLayout.Padding = [0 0 0 0];
      app.GridLayout.Scrollable = 'on';

      % Create LeftPanel
      app.LeftPanel = uipanel(app.GridLayout);
      app.LeftPanel.Layout.Row = 1;
      app.LeftPanel.Layout.Column = 1;

      % Create CADFileEditFieldLabel
      app.CADFileEditFieldLabel = uilabel(app.LeftPanel);
      app.CADFileEditFieldLabel.HorizontalAlignment = 'right';
      app.CADFileEditFieldLabel.Position = [8 341 54 22];
      app.CADFileEditFieldLabel.Text = 'CAD File';

      % Create CADFileEditField
      app.CADFileEditField = uieditfield(app.LeftPanel, 'text');
      app.CADFileEditField.Position = [77 341 126 22];

      % Create Button
      app.Button = uibutton(app.LeftPanel, 'push');
      app.Button.Position = [216 341 26 22];
      app.Button.Text = '...';

      % Create ScaleSpinnerLabel
      app.ScaleSpinnerLabel = uilabel(app.LeftPanel);
      app.ScaleSpinnerLabel.HorizontalAlignment = 'right';
      app.ScaleSpinnerLabel.Position = [94 269 35 22];
      app.ScaleSpinnerLabel.Text = 'Scale';

      % Create ScaleSpinner
      app.ScaleSpinner = uispinner(app.LeftPanel);
      app.ScaleSpinner.Position = [144 269 100 22];
      app.ScaleSpinner.Value = 1;
      app.ScaleSpinner.Step = 0.1;

      % Create XYMirrorSymmetryCheckBox
      app.XYMirrorSymmetryCheckBox = uicheckbox(app.LeftPanel);
      app.XYMirrorSymmetryCheckBox.Text = 'XY Mirror Symmetry';
      app.XYMirrorSymmetryCheckBox.Position = [11 103 130 22];

      % Create RotationalSymmetrySpinnerLabel
      app.RotationalSymmetrySpinnerLabel = uilabel(app.LeftPanel);
      app.RotationalSymmetrySpinnerLabel.HorizontalAlignment = 'right';
      app.RotationalSymmetrySpinnerLabel.Position = [8 136 117 22];
      app.RotationalSymmetrySpinnerLabel.Text = 'Rotational Symmetry';

      % Create RotationalSymmetrySpinner
      app.RotationalSymmetrySpinner = uispinner(app.LeftPanel);
      app.RotationalSymmetrySpinner.Limits = [0 Inf];
      app.RotationalSymmetrySpinner.Position = [140 136 63 22];
      app.RotationalSymmetrySpinner.Value = 1;

      % Create VariableNameEditFieldLabel
      app.VariableNameEditFieldLabel = uilabel(app.LeftPanel);
      app.VariableNameEditFieldLabel.HorizontalAlignment = 'right';
      app.VariableNameEditFieldLabel.Position = [8 379 84 22];
      app.VariableNameEditFieldLabel.Text = 'Variable Name';

      % Create VariableNameEditField
      app.VariableNameEditField = uieditfield(app.LeftPanel, 'text');
      app.VariableNameEditField.Position = [107 379 135 22];

      % Create ShowPreviewCheckBox
      app.ShowPreviewCheckBox = uicheckbox(app.LeftPanel);
      app.ShowPreviewCheckBox.Text = 'Show Preview';
      app.ShowPreviewCheckBox.Position = [8 43 98 22];
      app.ShowPreviewCheckBox.Value = true;

      % Create UpdateButton
      app.UpdateButton = uibutton(app.LeftPanel, 'push');
      app.UpdateButton.Enable = 'off';
      app.UpdateButton.Position = [158 14 84 22];
      app.UpdateButton.Text = 'Update';

      % Create OffsetSpinnerLabel
      app.OffsetSpinnerLabel = uilabel(app.LeftPanel);
      app.OffsetSpinnerLabel.HorizontalAlignment = 'right';
      app.OffsetSpinnerLabel.Position = [16 238 37 22];
      app.OffsetSpinnerLabel.Text = 'Offset';

      % Create OffsetSpinner
      app.OffsetSpinner = uispinner(app.LeftPanel);
      app.OffsetSpinner.Position = [60 238 65 22];

      % Create OffsetSpinner_2
      app.OffsetSpinner_2 = uispinner(app.LeftPanel);
      app.OffsetSpinner_2.Position = [124 238 59 22];

      % Create OffsetSpinner_3
      app.OffsetSpinner_3 = uispinner(app.LeftPanel);
      app.OffsetSpinner_3.Position = [182 238 60 22];

      % Create OffsetSpinner_4
      app.OffsetSpinner_4 = uispinner(app.LeftPanel);
      app.OffsetSpinner_4.Position = [124 206 59 22];

      % Create OffsetSpinner_5
      app.OffsetSpinner_5 = uispinner(app.LeftPanel);
      app.OffsetSpinner_5.Position = [182 206 60 22];

      % Create RotationSpinnerLabel
      app.RotationSpinnerLabel = uilabel(app.LeftPanel);
      app.RotationSpinnerLabel.HorizontalAlignment = 'right';
      app.RotationSpinnerLabel.Position = [3 206 50 22];
      app.RotationSpinnerLabel.Text = 'Rotation';

      % Create RotationSpinner
      app.RotationSpinner = uispinner(app.LeftPanel);
      app.RotationSpinner.Position = [60 206 65 22];

      % Create AutoupdateCheckBox
      app.AutoupdateCheckBox = uicheckbox(app.LeftPanel);
      app.AutoupdateCheckBox.Text = 'Auto-update';
      app.AutoupdateCheckBox.Position = [8 14 87 22];
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
      app.UIAxes.Position = [7 45 373 328];

      % Show the figure after all components are created
      app.UIFigure.Visible = 'on';
    end
  end

  methods (Access=public)
    function app = CadFileLoader()
      % Start the CadFileLoader interface

      % Create UI
      app.createComponents();
      
      % Register the app with App Designer
      registerApp(app, app.UIFigure);

      if nargout == 0
        clear app;
      end
    end

    % Code that executes before app deletion
    function delete(app)

      % Delete UIFigure when app is deleted
      delete(app.UIFigure)
    end
  end
end

