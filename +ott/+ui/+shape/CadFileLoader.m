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
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 420];
  end
  
  % Properties that correspond to app components
  properties (Access = public)
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
    OffsetXyzSpinners               ott.ui.support.XyzSpinners
    RotationXyzSpinners             ott.ui.support.XyzSpinners
    UIAxes                          matlab.ui.control.UIAxes
    updateCheckButton               ott.ui.support.UpdateCheckButton
  end

  % Component initialization
  methods (Access = protected)
    
    function startupFcn(app)
    end

    % Create UIFigure and components
    function createLeftComponents(app)
      
      lmargin = 10;

      % Create VariableNameEditFieldLabel
      app.VariableNameEditFieldLabel = uilabel(app.LeftPanel);
      app.VariableNameEditFieldLabel.HorizontalAlignment = 'left';
      app.VariableNameEditFieldLabel.Position = [lmargin 379 84 22];
      app.VariableNameEditFieldLabel.Text = 'Variable Name';

      % Create VariableNameEditField
      app.VariableNameEditField = uieditfield(app.LeftPanel, 'text');
      app.VariableNameEditField.Position = [107 379 135 22];

      % Create CADFileEditFieldLabel
      app.CADFileEditFieldLabel = uilabel(app.LeftPanel);
      app.CADFileEditFieldLabel.HorizontalAlignment = 'left';
      app.CADFileEditFieldLabel.Position = [lmargin 341 54 22];
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
      app.ScaleSpinnerLabel.HorizontalAlignment = 'left';
      app.ScaleSpinnerLabel.Position = [lmargin 269 35 22];
      app.ScaleSpinnerLabel.Text = 'Scale';

      % Create ScaleSpinner
      app.ScaleSpinner = uispinner(app.LeftPanel);
      app.ScaleSpinner.Position = [144 269 100 22];
      app.ScaleSpinner.Value = 1;
      app.ScaleSpinner.Step = 0.1;

      % Create XYMirrorSymmetryCheckBox
      app.XYMirrorSymmetryCheckBox = uicheckbox(app.LeftPanel);
      app.XYMirrorSymmetryCheckBox.Text = 'XY Mirror Symmetry';
      app.XYMirrorSymmetryCheckBox.Position = [lmargin 103 130 22];

      % Create RotationalSymmetrySpinnerLabel
      app.RotationalSymmetrySpinnerLabel = uilabel(app.LeftPanel);
      app.RotationalSymmetrySpinnerLabel.HorizontalAlignment = 'left';
      app.RotationalSymmetrySpinnerLabel.Position = [lmargin 136 117 22];
      app.RotationalSymmetrySpinnerLabel.Text = 'Rotational Symmetry';

      % Create RotationalSymmetrySpinner
      app.RotationalSymmetrySpinner = uispinner(app.LeftPanel);
      app.RotationalSymmetrySpinner.Limits = [0 Inf];
      app.RotationalSymmetrySpinner.Position = [140 136 63 22];
      app.RotationalSymmetrySpinner.Value = 1;

      % Create ShowPreviewCheckBox
      app.ShowPreviewCheckBox = uicheckbox(app.LeftPanel);
      app.ShowPreviewCheckBox.Text = 'Show Preview';
      app.ShowPreviewCheckBox.Position = [lmargin 43 98 22];
      app.ShowPreviewCheckBox.Value = true;
      
      % Offset spinners
      app.OffsetXyzSpinners = ott.ui.support.XyzSpinners(app.LeftPanel, ...
          'label', 'Offset', 'position', [10, 240]);
        
      % Rotation spinners
      app.RotationXyzSpinners = ott.ui.support.XyzSpinners(app.LeftPanel, ...
          'label', 'Rotation', 'position', [10, 210]);
      
      % Auto-update widget
      app.updateCheckButton = ott.ui.support.UpdateCheckButton(...
          app.LeftPanel, 'position', [lmargin, 14]);
        
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
      app.UIAxes.Position = [7 45 373 328];
      
    end
  end

  methods (Access=public)
    function app = CadFileLoader()
      % Start the CadFileLoader interface
      
      app = app@ott.ui.shape.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
end

