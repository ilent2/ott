classdef Simple < ott.ui.shape.NewShapeBase ...
    & ott.ui.support.GenerateCodeMenu
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

% TODO: It would be nice to change the pill tip height limits when
%   the pill height changes.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Simple Shape Generator';

    aboutText = ['Generates a simple geometric shape.  For more complex' ...
      ' shapes consider using a external CAD program and loading the' ...
      ' shape using CadFileLoader.'];
    
    helpText = {ott.ui.shape.Simple.aboutText, ...
      ''};
    
    windowName = ott.ui.shape.Simple.nameText;
    
    ShapeOptions = {'Sphere', 'Ellipsoid', 'Cylinder', ...
          'Cube', 'Rectangular Prism', 'Pill', 'Bicone', 'Cone Tipped Cylinder', ...
          'Biconcave Disc'};
    ShapeSmoothDraw = containers.Map(ott.ui.shape.Simple.ShapeOptions, ...
      [true, true, false, false, false, true, false, false, true]);
  end
  
  % Properties that correspond to app components
  properties (Access = public)
    ShapeDropDown               ott.ui.support.LabeledDropDown
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
  end
  
  methods (Access=private)
    function ShapeDropDownValueChanged(app, evt)
      
      % Update visible shape widgets
      app.UpdateVisibleShapeWidgets();
      
      % Continue processing changed value
      app.valueChangedCb(evt);
    end
      
    function UpdateVisibleShapeWidgets(app)
      % Update the visible widgets associated with the selected shape type
      
      rheight = 32;
      
      % Hide all widgets
      app.ExtraGrid.RowHeight(2:16) = repmat({0}, 1, 15);
      
      % Disable all widgets
      widgets = findobj(app.ExtraGrid.Children, '-property', 'Enable');
      set(widgets, 'Enable', 'off');
      app.ShapeDropDown.Enable = 'on';
      
      % Show widget for shape
      switch app.ShapeDropDown.Value
        case 'Sphere'
          app.ExtraGrid.RowHeight{app.SphereSpinner.Layout.Row} = rheight;
          app.SphereSpinner.Enable = 'on';
        case 'Ellipsoid'
          app.ExtraGrid.RowHeight{app.EllipsoidXyzSpinners.Layout.Row} = rheight;
          app.EllipsoidXyzSpinners.Enable = 'on';
        case 'Cylinder'
          app.ExtraGrid.RowHeight{app.CylinderRadiusSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.CylinderHeightSpinner.Layout.Row} = rheight;
          app.CylinderHeightSpinner.Enable = 'on';
          app.CylinderRadiusSpinner.Enable = 'on';
        case 'Cube'
          app.ExtraGrid.RowHeight{app.CubeSpinner.Layout.Row} = rheight;
          app.CubeSpinner.Enable = 'on';
        case 'Rectangular Prism'
          app.ExtraGrid.RowHeight{app.PrismXyzSpinners.Layout.Row} = rheight;
          app.PrismXyzSpinners.Enable = 'on';
        case 'Pill'
          app.ExtraGrid.RowHeight{app.PillRadiusSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.PillHeightSpinner.Layout.Row} = rheight;
          app.PillRadiusSpinner.Enable = 'on';
          app.PillHeightSpinner.Enable = 'on';
        case 'Bicone'
          app.ExtraGrid.RowHeight{app.BiconeRadiusSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.BiconeHeightSpinner.Layout.Row} = rheight;
          app.BiconeRadiusSpinner.Enable = 'on';
          app.BiconeHeightSpinner.Enable = 'on';
        case 'Cone Tipped Cylinder'
          app.ExtraGrid.RowHeight{app.ConeRadiusSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.ConeHeightSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.ConeTHeightSpinner.Layout.Row} = rheight;
          app.ConeRadiusSpinner.Enable = 'on';
          app.ConeHeightSpinner.Enable = 'on';
          app.ConeTHeightSpinner.Enable = 'on';
        case 'Biconcave Disc'
          app.ExtraGrid.RowHeight{app.DiscRadiusSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.DiscCoeffSpinners.Layout.Row} = rheight;
          app.DiscRadiusSpinner.Enable = 'on';
          app.DiscCoeffSpinners.Enable = 'on';
        otherwise
          error('Internal error');
      end
    end
  end
  
  methods (Access=protected)
    
    function setDefaultValues(app)
      
      scale = 1e-6;
      
      % App specifics
      app.ShapeDropDown.Value = 'Sphere';
      app.SphereSpinner.Value = 1*scale;
      app.EllipsoidXyzSpinners.Value = [1,1.2,1.4]*scale;
      app.CylinderRadiusSpinner.Value = 0.5*scale;
      app.CylinderHeightSpinner.Value = 1.0*scale;
      app.CubeSpinner.Value = 1.0*scale;
      app.PrismXyzSpinners.Value = [1,1.2,1.4]*scale;
      app.PillRadiusSpinner.Value = 0.5*scale;
      app.PillHeightSpinner.Value = 2.0*scale;
      app.BiconeRadiusSpinner.Value = 0.5*scale;
      app.BiconeHeightSpinner.Value = 1.0*scale;
      app.ConeRadiusSpinner.Value = 0.5*scale;
      app.ConeHeightSpinner.Value = 1.0*scale;
      app.ConeTHeightSpinner.Value = 0.25*scale;
      app.DiscRadiusSpinner.Value = 7.82*scale;
      app.DiscCoeffSpinners.Value = [0.0518, 2.0026, -4.491];
      
      % Base class/finalisation
      % TODO: This is kind of kludgy, we update the window twice!!!
      setDefaultValues@ott.ui.shape.NewShapeBase(app);
      
      % Display appropriate widgets
      app.UpdateVisibleShapeWidgets();
    end
    
    function data = GenerateData(app)
      % Generate shape
      
      switch app.ShapeDropDown.Value
        case 'Sphere'
          data = ott.shape.Sphere(app.SphereSpinner.Value);
        case 'Ellipsoid'
          data = ott.shape.Ellipsoid(app.EllipsoidXyzSpinners.Value);
        case 'Cylinder'
          data = ott.shape.Cylinder(app.CylinderRadiusSpinner.Value, ...
            app.CylinderHeightSpinner.Value);
        case 'Cube'
          data = ott.shape.Cube(app.CubeSpinner.Value);
        case 'Rectangular Prism'
          data = ott.shape.RectangularPrism(app.PrismXyzSpinners.Value);
        case 'Pill'
          data = ott.shape.AxisymFunc.Pill(app.PillHeightSpinner.Value, ...
            app.PillRadiusSpinner.Value);
        case 'Bicone'
          data = ott.shape.AxisymInterp.Bicone(app.BiconeHeightSpinner.Value, ...
            app.BiconeRadiusSpinner.Value);
        case 'Cone Tipped Cylinder'
          data = ott.shape.AxisymInterp.ConeTippedCylinder(...
            app.ConeHeightSpinner.Value, app.ConeRadiusSpinner.Value, ...
            app.ConeTHeightSpinner.Value);
        case 'Biconcave Disc'
          data = ott.shape.AxisymFunc.BiconcaveDisc(...
            app.DiscRadiusSpinner.Value, app.DiscCoeffSpinners.Value);
        otherwise
          error('Internal error');
      end
    end
    
    function code = generateCode(app)
      % Generate visualisation code
      
      code = {};
      
      % Shape code
      code{end+1} = '% Generate shape';
      switch app.ShapeDropDown.Value
        case 'Sphere'
          code{end+1} = ['radius = ' num2str(app.SphereSpinner.Value) ';'];
          code{end+1} = 'shape = ott.shape.Sphere(radius);';
        case 'Ellipsoid'
          code{end+1} = ['xyzSize = ' num2str(app.EllipsoidXyzSpinners.Value) ';'];
          code{end+1} = 'shape = ott.shape.Ellipsoid(xyzSize);';
        case 'Cylinder'
          code{end+1} = ['radius = ' num2str(app.CylinderRadiusSpinner.Value) ';'];
          code{end+1} = ['height = ' num2str(app.CylinderHeightSpinner.Value) ';'];
          code{end+1} = 'shape = ott.shape.Cylinder(radius, height);';
        case 'Cube'
          code{end+1} = ['width = ' num2str(app.CubeSpinner.Value) ';'];
          code{end+1} = 'shape = ott.shape.Cube(width);';
        case 'Rectangular Prism'
          code{end+1} = ['xyzSize = [' num2str(app.PrismXyzSpinners.Value) '];'];
          code{end+1} = 'shape = ott.shape.RectangularPrism(xyzSize);';
        case 'Pill'
          code{end+1} = ['height = ' num2str(app.PillHeightSpinner.Value) ';'];
          code{end+1} = ['radius = ' num2str(app.PillRadiusSpinner.Value) ';'];
          code{end+1} = 'shape = ott.shape.AxisymFunc.Pill(height, radius);';
        case 'Bicone'
          code{end+1} = ['height = ' num2str(app.BiconeHeightSpinner.Value) ';'];
          code{end+1} = ['radius = ' num2str(app.BiconeRadiusSpinner.Value) ';'];
          code{end+1} = 'shape = ott.shape.AxisymInterp.Bicone(height, radius);';
        case 'Cone Tipped Cylinder'
          code{end+1} = ['height = ' num2str(app.ConeHeightSpinner.Value) ';'];
          code{end+1} = ['radius = ' num2str(app.ConeRadiusSpinner.Value) ';'];
          code{end+1} = ['tipHeight = ' num2str(app.ConeTHeightSpinner.Value) ';'];
          code{end+1} = 'shape = ott.shape.AxisymInterp.ConeTippedCylinder(...';
          code{end+1} = '    height, radius, tipHeight);';
        case 'Biconcave Disc'
          code{end+1} = ['radius = ' num2str(app.DiscRadiusSpinner.Value) ';'];
          code{end+1} = ['coefficients = [' num2str(app.DiscCoeffSpinners.Value) '];'];
          code{end+1} = 'shape = ott.shape.AxisymFunc.BiconcaveDisc(radius, coefficients);';
        otherwise
          error('Internal error');
      end
      
      % Clear line
      code{end+1} = '';
      code{end+1} = '% Generate visualisation';
      
      % Visualisation code
      code{end+1} = 'figure()';
      
      if app.ShapeSmoothDraw(app.ShapeDropDown.Value)
        code{end+1} = 'shape.surf(''surfOptions'', {''EdgeColor'', ''none''});';
        code{end+1} = 'camlight; lighting gouraud;';
      else
        code{end+1} = 'shape.surf();';
      end
    end
    
    function updateShapePreview(app)
      
      % Call base update
      updateShapePreview@ott.ui.shape.NewShapeBase(app);
      
      % Change visualisation style for specific shapes
      if app.ShapeSmoothDraw(app.ShapeDropDown.Value)
        camlight(app.PreviewUIAxes);
        lighting(app.PreviewUIAxes, 'gouraud');
        pt = findobj(app.PreviewUIAxes.Children, '-class', 'matlab.graphics.primitive.Patch');
        pt.EdgeColor = 'none';
      end
      
    end
    
    function createLeftComponents(app)
      % Create additional components for application
      
      % Call base for most things
      createLeftComponents@ott.ui.shape.NewShapeBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 17);
      app.ExtraGrid.RowHeight{end} = '1x';
      
      % Create ShapeDropDown
      app.ShapeDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
          'label', 'Shape');
      app.ShapeDropDown.Items = app.ShapeOptions;
      app.ShapeDropDown.ValueChangedFcn = ...
          createCallbackFcn(app, @ShapeDropDownValueChanged, true);
      app.ShapeDropDown.Value = 'Sphere';
      app.ShapeDropDown.Layout.Row = 1;
      app.ShapeDropDown.Layout.Column = 1;
      
      % Create components for each shape, make hidden by default
      app.SphereSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Radius');
      app.SphereSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.SphereSpinner.Step = 1e-7;
      app.SphereSpinner.Limits = [0, Inf];
      app.SphereSpinner.LowerLimitInclusive = 'off';
      app.SphereSpinner.Layout.Row = 2;
      app.SphereSpinner.Layout.Column = 1;
      
      app.EllipsoidXyzSpinners = ott.ui.support.XyzSpinners(app.ExtraGrid, ...
          'label', 'Radii');
      app.EllipsoidXyzSpinners.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.EllipsoidXyzSpinners.Step = 1e-7;
      app.EllipsoidXyzSpinners.Limits = [0, Inf];
      app.EllipsoidXyzSpinners.LowerLimitInclusive = 'off';
      app.EllipsoidXyzSpinners.Layout.Row = 3;
      app.EllipsoidXyzSpinners.Layout.Column = 1;
      
      app.CylinderRadiusSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Radius');
      app.CylinderRadiusSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.CylinderRadiusSpinner.Step = 1e-7;
      app.CylinderRadiusSpinner.Limits = [0, Inf];
      app.CylinderRadiusSpinner.LowerLimitInclusive = 'off';
      app.CylinderRadiusSpinner.Layout.Row = 4;
      app.CylinderRadiusSpinner.Layout.Column = 1;
      
      app.CylinderHeightSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Height');
      app.CylinderHeightSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.CylinderHeightSpinner.Step = 1e-7;
      app.CylinderHeightSpinner.Limits = [0, Inf];
      app.CylinderHeightSpinner.LowerLimitInclusive = 'off';
      app.CylinderHeightSpinner.Layout.Row = 5;
      app.CylinderHeightSpinner.Layout.Column = 1;
      
      app.CubeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Width');
      app.CubeSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.CubeSpinner.Step = 1e-7;
      app.CubeSpinner.Limits = [0, Inf];
      app.CubeSpinner.LowerLimitInclusive = 'off';
      app.CubeSpinner.Layout.Row = 6;
      app.CubeSpinner.Layout.Column = 1;
      
      app.PrismXyzSpinners = ott.ui.support.XyzSpinners(app.ExtraGrid, ...
          'label', 'Width');
      app.PrismXyzSpinners.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.PrismXyzSpinners.Step = 1e-7;
      app.PrismXyzSpinners.Limits = [0, Inf];
      app.PrismXyzSpinners.LowerLimitInclusive = 'off';
      app.PrismXyzSpinners.Layout.Row = 7;
      app.PrismXyzSpinners.Layout.Column = 1;
      
      app.PillRadiusSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Radius');
      app.PillRadiusSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.PillRadiusSpinner.Step = 1e-7;
      app.PillRadiusSpinner.Limits = [0, Inf];
      app.PillRadiusSpinner.LowerLimitInclusive = 'off';
      app.PillRadiusSpinner.Layout.Row = 8;
      app.PillRadiusSpinner.Layout.Column = 1;
      
      app.PillHeightSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Height');
      app.PillHeightSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.PillHeightSpinner.Step = 1e-7;
      app.PillHeightSpinner.Limits = [0, Inf];
      app.PillHeightSpinner.LowerLimitInclusive = 'off';
      app.PillHeightSpinner.Layout.Row = 9;
      app.PillHeightSpinner.Layout.Column = 1;
      
      app.BiconeRadiusSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Radius');
      app.BiconeRadiusSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.BiconeRadiusSpinner.Step = 1e-7;
      app.BiconeRadiusSpinner.Limits = [0, Inf];
      app.BiconeRadiusSpinner.LowerLimitInclusive = 'off';
      app.BiconeRadiusSpinner.Layout.Row = 10;
      app.BiconeRadiusSpinner.Layout.Column = 1;
      
      app.BiconeHeightSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Height');
      app.BiconeHeightSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.BiconeHeightSpinner.Step = 1e-7;
      app.BiconeHeightSpinner.Limits = [0, Inf];
      app.BiconeHeightSpinner.LowerLimitInclusive = 'off';
      app.BiconeHeightSpinner.Layout.Row = 11;
      app.BiconeHeightSpinner.Layout.Column = 1;
      
      app.ConeRadiusSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Radius');
      app.ConeRadiusSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.ConeRadiusSpinner.Step = 1e-7;
      app.ConeRadiusSpinner.Limits = [0, Inf];
      app.ConeRadiusSpinner.LowerLimitInclusive = 'off';
      app.ConeRadiusSpinner.Layout.Row = 12;
      app.ConeRadiusSpinner.Layout.Column = 1;
      
      app.ConeHeightSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Total Height');
      app.ConeHeightSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.ConeHeightSpinner.Step = 1e-7;
      app.ConeHeightSpinner.Limits = [0, Inf];
      app.ConeHeightSpinner.LowerLimitInclusive = 'off';
      app.ConeHeightSpinner.Layout.Row = 13;
      app.ConeHeightSpinner.Layout.Column = 1;
      
      app.ConeTHeightSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Tip Height');
      app.ConeTHeightSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.ConeTHeightSpinner.Step = 1e-7;
      app.ConeTHeightSpinner.Limits = [0, Inf];
      app.ConeTHeightSpinner.LowerLimitInclusive = 'off';
      app.ConeTHeightSpinner.Layout.Row = 14;
      app.ConeTHeightSpinner.Layout.Column = 1;
      
      app.DiscRadiusSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Radius');
      app.DiscRadiusSpinner.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.DiscRadiusSpinner.Step = 1e-7;
      app.DiscRadiusSpinner.Limits = [0, Inf];
      app.DiscRadiusSpinner.LowerLimitInclusive = 'off';
      app.DiscRadiusSpinner.Layout.Row = 15;
      app.DiscRadiusSpinner.Layout.Column = 1;
      
      app.DiscCoeffSpinners = ott.ui.support.XyzSpinners(app.ExtraGrid, ...
          'label', 'Coefficients');
      app.DiscCoeffSpinners.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      app.DiscCoeffSpinners.Step = 0.1;
      app.DiscCoeffSpinners.Layout.Row = 16;
      app.DiscCoeffSpinners.Layout.Column = 1;
      
    end
  end
  
  methods (Access=public)
    function app=Simple()
      % Start the simple shape creation GUI
      
      app = app@ott.ui.shape.NewShapeBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end