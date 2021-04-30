classdef Simple < ott.ui.support.AppTopLevel
% Constructs a drag tensor for a shape.
% This interface uses :meth:`ott.drag.Stokes.FromShape`.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Simple Drag';

    aboutText = ['Constructs a drag tensor for a shape.  This interface' ...
      'uses ott.drag.Stokes.FromShape.  Not all shapes are supported, ' ...
      'most shapes default to a sphere whose radius matches the ' ...
      'maxRadius of the shape.'];
    
    helpText = {ott.ui.drag.Simple.aboutText, ...
      ''};
    
    windowName = ott.ui.drag.Simple.nameText;
    windowSize = [550, 200];
  end
  
  % Properties that correspond to app components
  properties (Access = public)
    VariableNameEditFieldLabel          matlab.ui.control.Label
    VariableNameEditField               matlab.ui.control.EditField
    ShapeDropDownLabel                  matlab.ui.control.Label
    ShapeDropDown                       matlab.ui.control.DropDown
    ViscositySpinnerLabel               matlab.ui.control.Label
    ViscositySpinner                    matlab.ui.control.Spinner
    UITable                             matlab.ui.control.Table
    updateCheckButton               ott.ui.support.UpdateCheckButton
  end

  methods (Access=protected)
    function startupFcn(app)
    end
    
    function createMainComponents(app)
      
      lmargin = 10;
      wsize = [135 22];
      x0 = 107;

      % Create VariableNameEditFieldLabel
      app.VariableNameEditFieldLabel = uilabel(app.UIFigure);
      app.VariableNameEditFieldLabel.HorizontalAlignment = 'left';
      app.VariableNameEditFieldLabel.Position = [lmargin 158 84 22];
      app.VariableNameEditFieldLabel.Text = 'Variable Name';

      % Create VariableNameEditField
      app.VariableNameEditField = uieditfield(app.UIFigure, 'text');
      app.VariableNameEditField.Position = [x0 158 wsize];

      % Create ShapeDropDownLabel
      app.ShapeDropDownLabel = uilabel(app.UIFigure);
      app.ShapeDropDownLabel.HorizontalAlignment = 'left';
      app.ShapeDropDownLabel.Position = [lmargin 128 40 22];
      app.ShapeDropDownLabel.Text = 'Shape';

      % Create ShapeDropDown
      app.ShapeDropDown = uidropdown(app.UIFigure);
      app.ShapeDropDown.Position = [x0 128 wsize];

      % Create ViscositySpinnerLabel
      app.ViscositySpinnerLabel = uilabel(app.UIFigure);
      app.ViscositySpinnerLabel.HorizontalAlignment = 'left';
      app.ViscositySpinnerLabel.Position = [lmargin 98 52 22];
      app.ViscositySpinnerLabel.Text = 'Viscosity';

      % Create ViscositySpinner
      app.ViscositySpinner = uispinner(app.UIFigure);
      app.ViscositySpinner.Position = [x0 98 wsize];
      
      % Create UITable
      cwidth = 40;
      app.UITable = uitable(app.UIFigure);
      app.UITable.ColumnName = {'x', 'y', 'z', 'Rx', 'Ry', 'Rz'};
      app.UITable.RowName = {'x', 'y', 'z', 'Rx', 'Ry', 'Rz'};
      app.UITable.ColumnWidth = repmat({cwidth}, 1, 6);
      app.UITable.Data = repmat({''}, 6, 6);
      app.UITable.Position = [250 20 cwidth*6+42 160];
      
      % Auto-update widget
      app.updateCheckButton = ott.ui.support.UpdateCheckButton(...
          app.UIFigure, 'position', [lmargin, 20]);
    end
  end
  
  methods (Access=public)
    function app=Simple()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTopLevel();
      if nargout == 0
        clear app;
      end
    end
  end
  
end