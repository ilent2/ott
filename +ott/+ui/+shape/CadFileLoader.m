classdef CadFileLoader < ott.ui.shape.NewShapeBase ...
    & ott.ui.support.GenerateCodeMenu
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
    
    windowName = ott.ui.shape.CadFileLoader.nameText;
  end
  
  % Properties that correspond to app components
  properties (Access = public)
    CadFileSelector                 ott.ui.support.FileSelector
    ScaleSpinner                    ott.ui.support.LabeledSpinner
    XYMirrorSymmetryCheckBox        matlab.ui.control.CheckBox
    RotationalSymmetrySpinner       ott.ui.support.LabeledSpinner
  end

  % Component initialization
  methods (Access = protected)
    
    function setDefaultValues(app)
      
      % App specifics
      app.CadFileSelector.Value = '';
      app.ScaleSpinner.Value = 1;
      app.XYMirrorSymmetryCheckBox.Value = true;
      app.RotationalSymmetrySpinner.Value = 1;
      
      % Base class/finalisation
      setDefaultValues@ott.ui.shape.NewShapeBase(app);
    end
    
    function shape = generateShape(app)
      % TODO
      shape = [];
    end
    
    function code = generateCode(app)
      % TODO
      code = {};
    end

    % Create UIFigure and components
    function createLeftComponents(app)
      % Create additional components for application
      
      % Call base for most things
      createLeftComponents@ott.ui.shape.NewShapeBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = {32, 32, 32, 32, '1x'};
      
      % Create CAD file selector
      app.CadFileSelector = ott.ui.support.FileSelector(app.ExtraGrid, ...
        'label', 'CAD File');
      app.CadFileSelector.Filter = {'*.obj;*.stl', 'CAD file'};
      app.CadFileSelector.Title = 'Select a File';
      app.CadFileSelector.PostSelect = app.UIFigure;
      app.CadFileSelector.Layout.Row = 1;
      app.CadFileSelector.Layout.Column = 1;

      % Create ScaleSpinnerLabel
      app.ScaleSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Scale');
      app.ScaleSpinner.Step = 0.1;
      app.ScaleSpinner.Layout.Row = 2;
      app.ScaleSpinner.Layout.Column = 1;

      % Create XYMirrorSymmetryCheckBox
      app.XYMirrorSymmetryCheckBox = uicheckbox(app.ExtraGrid);
      app.XYMirrorSymmetryCheckBox.Text = 'XY Mirror Symmetry';
      app.XYMirrorSymmetryCheckBox.Layout.Row = 3;
      app.XYMirrorSymmetryCheckBox.Layout.Column = 1;

      % Create RotationalSymmetrySpinnerLabel
      app.RotationalSymmetrySpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Rotational Symmetry');
      app.RotationalSymmetrySpinner.Limits = [0 Inf];
      app.RotationalSymmetrySpinner.Layout.Row = 4;
      app.RotationalSymmetrySpinner.Layout.Column = 1;
        
    end
  end

  methods (Access=public)
    function app = CadFileLoader()
      % Start the CadFileLoader interface
      
      app = app@ott.ui.shape.NewShapeBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
end

