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
%   shape = app.Data

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
      ['Rotational/XY Symmetry/starShaped: Set these if the objects has', ...
      'symetries (including rotations/translations).  If in ', ...
      'doubt, do not set these.']};
    
    windowName = ott.ui.shape.CadFileLoader.nameText;
  end
  
  % Properties that correspond to app components
  properties (Access = public)
    CadFileSelector                 ott.ui.support.FileSelector
    ScaleSpinner                    ott.ui.support.LabeledSpinner
    XYMirrorSymmetryCheckBox        matlab.ui.control.CheckBox
    StarShapedCheckBox              matlab.ui.control.CheckBox
    RotationalSymmetrySpinner       ott.ui.support.LabeledSpinner
  end

  % Component initialization
  methods (Access = protected)
    
    function setDefaultValues(app)
      
      % App specifics
      app.CadFileSelector.Value = '';
      app.ScaleSpinner.Value = 1;
      app.XYMirrorSymmetryCheckBox.Value = false;
      app.StarShapedCheckBox.Value = false;
      app.RotationalSymmetrySpinner.Value = 1;
      
      % Base class/finalisation
      setDefaultValues@ott.ui.shape.NewShapeBase(app, false);
    end
    
    function data = generateData(app)
      % Load a CAD file
      
      fname = app.CadFileSelector.Value;
      if isempty(fname)
        data = [];
        warning('ott:ui:shape:CadFileLoader:no_filename', 'No shape specified');
        return;
      end
      
      [~, ~, ext] = fileparts(fname);
      
      % Get scale
      scale = app.ScaleSpinner.Value;
      
      % Get rotation/offset
      offset = app.OffsetXyzSpinners.Value(:) ./ scale;
      rotation = ott.utils.euler2rot(app.RotationXyzSpinners.Value(:));
      
      switch ext
        case '.obj'
          data = ott.shape.ObjLoader(fname, ...
            'position', offset, 'rotation', rotation);
        case '.stl'
          data = ott.shape.StlLoader(fname, ...
            'position', offset, 'rotation', rotation);
        otherwise
          error('Unsupported CAD file, must be .obj or .stl');
      end
      
      % Apply scale
      data = data * scale;
      
      % Set symmetry parameters
      data.starShaped = app.StarShapedCheckBox.Value;
      data.xySymmetry = app.XYMirrorSymmetryCheckBox.Value;
      data.zRotSymmetry = app.RotationalSymmetrySpinner.Value;
    end
    
    function code = generateCode(app)
      
      fname = app.CadFileSelector.Value;
      if isempty(fname)
        error('Must specify file name to generate code');
      end
      
      code = {};
      code{end+1} = ['fname = ''', fname, ''';'];
      code{end+1} = '';
      
      [~, ~, ext] = fileparts(fname);
      
      % Get scale
      code{end+1} = ['scale = ', num2str(app.ScaleSpinner.Value), ';'];
      
      % Get rotation/offset
      code{end+1} = ['offset = [', num2str(app.OffsetXyzSpinners.Value), '].'' ./ scale;'];
      code{end+1} = ['rotation = ott.utils.euler2rot([', num2str(app.RotationXyzSpinners.Value), '].'');'];
      code{end+1} = '';
      
      switch ext
        case '.obj'
          code{end+1} = 'shape = ott.shape.ObjLoader(fname, ...';
          code{end+1} = '  ''position'', offset, ''rotation'', rotation);';
        case '.stl'
          code{end+1} = 'shape = ott.shape.StlLoader(fname, ...';
          code{end+1} = '  ''position'', offset, ''rotation'', rotation);';
        otherwise
          error('Unsupported CAD file, must be .obj or .stl');
      end
      code{end+1} = 'shape = shape * scale;';
      code{end+1} = '';
      
      if app.StarShapedCheckBox.Value
        code{end+1} = 'shape.starShaped = true;';
      else
        code{end+1} = 'shape.starShaped = false;';
      end
      
      if app.XYMirrorSymmetryCheckBox.Value
        code{end+1} = 'shape.xySymmetry = true;';
      else
        code{end+1} = 'shape.xySymmetry = false;';
      end
      
      code{end+1} = ['shape.zRotSymmetry = ' num2str(app.RotationalSymmetrySpinner.Value), ';'];
    end

    % Create UIFigure and components
    function createLeftComponents(app)
      % Create additional components for application
      
      % Call base for most things
      createLeftComponents@ott.ui.shape.NewShapeBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 6);
      app.ExtraGrid.RowHeight{end} = '1x';
      
      % Create CAD file selector
      app.CadFileSelector = ott.ui.support.FileSelector(app.ExtraGrid, ...
        'label', 'CAD File');
      app.CadFileSelector.Filter = {'*.obj;*.stl', 'CAD file'};
      app.CadFileSelector.Title = 'Select a File';
      app.CadFileSelector.PostSelect = app.UIFigure;
      app.CadFileSelector.Layout.Row = 1;
      app.CadFileSelector.Layout.Column = 1;
      app.CadFileSelector.ValueChangedFcn = @(~,~) app.updateParametersCb();

      % Create ScaleSpinnerLabel
      app.ScaleSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Scale');
      app.ScaleSpinner.Step = 0.1;
      app.ScaleSpinner.Layout.Row = 2;
      app.ScaleSpinner.Layout.Column = 1;
      app.ScaleSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();

      % Create XYMirrorSymmetryCheckBox
      app.XYMirrorSymmetryCheckBox = uicheckbox(app.ExtraGrid);
      app.XYMirrorSymmetryCheckBox.Text = 'XY Mirror Symmetry';
      app.XYMirrorSymmetryCheckBox.Layout.Row = 3;
      app.XYMirrorSymmetryCheckBox.Layout.Column = 1;
      app.XYMirrorSymmetryCheckBox.ValueChangedFcn = @(~,~) app.updateParametersCb();

      % Create StarShapedCheckBox
      app.StarShapedCheckBox = uicheckbox(app.ExtraGrid);
      app.StarShapedCheckBox.Text = 'Star Shaped';
      app.StarShapedCheckBox.Layout.Row = 4;
      app.StarShapedCheckBox.Layout.Column = 1;
      app.StarShapedCheckBox.ValueChangedFcn = @(~,~) app.updateParametersCb();

      % Create RotationalSymmetrySpinner
      app.RotationalSymmetrySpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Rotational Symmetry');
      app.RotationalSymmetrySpinner.Limits = [0 Inf];
      app.RotationalSymmetrySpinner.Layout.Row = 5;
      app.RotationalSymmetrySpinner.Layout.Column = 1;
      app.RotationalSymmetrySpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
        
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

