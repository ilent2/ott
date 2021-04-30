classdef Simple < ott.ui.shape.AppBase
% Generate a simple geometric shape.
% For more complex shapes, consider using an external CAD program
% and :class:`CadFileLoader`.
%
% Supported shapes (and paramters):
%   - Sphere (radius)
%   - Ellipsoid (xyz-radius)
%   - Cylinder (radius, height)
%   - Cube (width)
%   - Rectangular Prism (xyz-width)
%   - Pill (radius, height)
%   - Bicone (radius, height)
%   - Cone Tipped Cylinder (radius, height, theight)
%   - Biconcave Disc (radius, 3-coefficients)
%
% Some of these shapes may move to their own GUI in a future release.
%
% This GUI can be launched from the launcher under
% Shape -> Simple or running the following command:
%
%   ott.ui.shape.Simple()
%
% The shape is stored internally and/or written to the matlab workspace
% if a variable name is given for the shape.  To access the internal
% instance use:
%
%   app = ott.ui.shape.Simple()
%   shape = app.shape

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Simple Shape Generator';

    aboutText = ['Generates a simple geometric shape.  For more complex' ...
      ' shapes consider using a external CAD program and loading the' ...
      ' shape using CadFileLoader.'];
    
    helpText = {ott.ui.shape.Simple.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 370];
  end
  
  % Properties that correspond to app components
  properties (Access = public)
    ShapeDropDownLabel          matlab.ui.control.Label
    ShapeDropDown               matlab.ui.control.DropDown
    VariableNameEditFieldLabel  matlab.ui.control.Label
    VariableNameEditField       matlab.ui.control.EditField
    SphereSpinner               ott.ui.support.LabeledSpinner
    EllipsoidXyzSpinners        ott.ui.support.XyzSpinners
    CylinderRadiusSpinner       ott.ui.support.LabeledSpinner
    CylinderHeightSpinner       ott.ui.support.LabeledSpinner
    CubeSpinner                 ott.ui.support.LabeledSpinner
    PrismXyzSpinners            ott.ui.support.XyzSpinners
    PillRadiusSpinner           ott.ui.support.LabeledSpinner
    PillHeightSpinner           ott.ui.support.LabeledSpinner
    BiconeRadiusSpinner         ott.ui.support.LabeledSpinner
    BiconeHeightSpinner         ott.ui.support.LabeledSpinner
    ConeRadiusSpinner           ott.ui.support.LabeledSpinner
    ConeHeightSpinner           ott.ui.support.LabeledSpinner
    ConeTHeightSpinner          ott.ui.support.LabeledSpinner
    DiscRadiusSpinner           ott.ui.support.LabeledSpinner
    DiscCoeffSpinners           ott.ui.support.XyzSpinners
    updateCheckButton           ott.ui.support.UpdateCheckButton
    ShowPreviewCheckBox         matlab.ui.control.CheckBox
  end
  
  methods (Access=private)
    function ShapeDropDownValueChanged(app, ~)
      
      % Hide all widgets
      app.SphereSpinner.Visible = 'off';
      app.EllipsoidXyzSpinners.Visible = 'off';
      app.CylinderRadiusSpinner.Visible = 'off';
      app.CylinderHeightSpinner.Visible = 'off';
      app.CubeSpinner.Visible = 'off';
      app.PrismXyzSpinners.Visible = 'off';
      app.PillRadiusSpinner.Visible = 'off';
      app.PillHeightSpinner.Visible = 'off';
      app.BiconeRadiusSpinner.Visible = 'off';
      app.BiconeHeightSpinner.Visible = 'off';
      app.ConeRadiusSpinner.Visible = 'off';
      app.ConeHeightSpinner.Visible = 'off';
      app.ConeTHeightSpinner.Visible = 'off';
      app.DiscRadiusSpinner.Visible = 'off';
      app.DiscCoeffSpinners.Visible = 'off';
      
      % Show widget for shape
      switch app.ShapeDropDown.Value
        case 'Sphere'
          app.SphereSpinner.Visible = 'on';
        case 'Ellipsoid'
          app.EllipsoidXyzSpinners.Visible = 'on';
        case 'Cylinder'
          app.CylinderRadiusSpinner.Visible = 'on';
          app.CylinderHeightSpinner.Visible = 'on';
        case 'Cube'
          app.CubeSpinner.Visible = 'on';
        case 'Rectangular Prism'
          app.PrismXyzSpinners.Visible = 'on';
        case 'Pill'
          app.PillRadiusSpinner.Visible = 'on';
          app.PillHeightSpinner.Visible = 'on';
        case 'Bicone'
          app.BiconeRadiusSpinner.Visible = 'on';
          app.BiconeHeightSpinner.Visible = 'on';
        case 'Cone Tipped Cylinder'
          app.ConeRadiusSpinner.Visible = 'on';
          app.ConeHeightSpinner.Visible = 'on';
          app.ConeTHeightSpinner.Visible = 'on';
        case 'Biconcave Disc'
          app.DiscRadiusSpinner.Visible = 'on';
          app.DiscCoeffSpinners.Visible = 'on';
        otherwise
          error('Internal error');
      end
    end
  end
  
  methods (Access=protected)
    function startupFcn(app)
      
      % Display appropriate widgets
      app.ShapeDropDownValueChanged();
      
    end
    
    function createLeftComponents(app)
      
      lmargin = 10;

      % Create VariableNameEditFieldLabel
      app.VariableNameEditFieldLabel = uilabel(app.LeftPanel);
      app.VariableNameEditFieldLabel.HorizontalAlignment = 'left';
      app.VariableNameEditFieldLabel.Position = [lmargin 330 84 22];
      app.VariableNameEditFieldLabel.Text = 'Variable Name';

      % Create VariableNameEditField
      app.VariableNameEditField = uieditfield(app.LeftPanel, 'text');
      app.VariableNameEditField.Position = [110 330 130 22];
      
      % Create ShapeDropDownLabel
      app.ShapeDropDownLabel = uilabel(app.LeftPanel);
      app.ShapeDropDownLabel.HorizontalAlignment = 'left';
      app.ShapeDropDownLabel.Position = [lmargin 300 40 22];
      app.ShapeDropDownLabel.Text = 'Shape';

      % Create ShapeDropDown
      app.ShapeDropDown = uidropdown(app.LeftPanel);
      app.ShapeDropDown.Items = {'Sphere', 'Ellipsoid', 'Cylinder', ...
          'Cube', 'Rectangular Prism', 'Pill', 'Bicone', 'Cone Tipped Cylinder', ...
          'Biconcave Disc'};
      app.ShapeDropDown.ValueChangedFcn = ...
          createCallbackFcn(app, @ShapeDropDownValueChanged, true);
      app.ShapeDropDown.Position = [110 300 130 22];
      app.ShapeDropDown.Value = 'Sphere';
      
      y1 = 250;
      y2 = y1 - 30;
      y3 = y2 - 30;
      
      % Create components for each shape, make hidden by default
      app.SphereSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Radius', 'position', [lmargin, y1], 'visible', 'off');
      app.EllipsoidXyzSpinners = ott.ui.support.XyzSpinners(app.LeftPanel, ...
          'label', 'Radii', 'position', [lmargin, y1], 'visible', 'off');
      app.CylinderRadiusSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Radius', 'position', [lmargin, y1], 'visible', 'off');
      app.CylinderHeightSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Height', 'position', [lmargin, y2], 'visible', 'off');
      app.CubeSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Width', 'position', [lmargin, y1], 'visible', 'off');
      app.PrismXyzSpinners = ott.ui.support.XyzSpinners(app.LeftPanel, ...
          'label', 'Width', 'position', [lmargin, y1], 'visible', 'off');
      app.PillRadiusSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Radius', 'position', [lmargin, y1], 'visible', 'off');
      app.PillHeightSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Height', 'position', [lmargin, y2], 'visible', 'off');
      app.BiconeRadiusSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Radius', 'position', [lmargin, y1], 'visible', 'off');
      app.BiconeHeightSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Height', 'position', [lmargin, y2], 'visible', 'off');
      app.ConeRadiusSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Radius', 'position', [lmargin, y1], 'visible', 'off');
      app.ConeHeightSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Height', 'position', [lmargin, y2], 'visible', 'off');
      app.ConeTHeightSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Tip Height', 'position', [lmargin, y3], 'visible', 'off');
      app.DiscRadiusSpinner = ott.ui.support.LabeledSpinner(app.LeftPanel, ...
          'label', 'Radius', 'position', [lmargin, y1], 'visible', 'off');
      app.DiscCoeffSpinners = ott.ui.support.XyzSpinners(app.LeftPanel, ...
          'label', 'Coefficients', 'position', [lmargin, y2], 'visible', 'off');
      
      % Create ShowPreviewCheckBox
      app.ShowPreviewCheckBox = uicheckbox(app.LeftPanel);
      app.ShowPreviewCheckBox.Text = 'Show Preview';
      app.ShowPreviewCheckBox.Position = [lmargin 43 98 22];
      app.ShowPreviewCheckBox.Value = true;
      
      % Auto-update widget
      app.updateCheckButton = ott.ui.support.UpdateCheckButton(...
          app.LeftPanel, 'position', [lmargin, 14]);
    end
    
    function createRightComponents(app)
    end
  end
  
  methods (Access=public)
    function app=Simple()
      % Start the simple shape creation GUI
      
      app = app@ott.ui.shape.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end