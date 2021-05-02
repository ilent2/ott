classdef Simple < ott.ui.support.AppTopLevel ...
    & ott.ui.support.RefreshInputsMenu ...
    & ott.ui.support.GenerateCodeMenu
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
    windowSize = [522, 179];
  end
  
  % Properties that correspond to app components
  properties (Access = public)
    
    MainGrid              matlab.ui.container.GridLayout
    
    % Left column
    VariableName          ott.ui.support.OutputVariableEntry
    ShapeDropdown         ott.ui.support.VariableDropdown
    ViscositySpinner      ott.ui.support.LabeledSpinner
    updateCheckButton     ott.ui.support.UpdateCheckButton
    
    % Right column
    UITable                             matlab.ui.control.Table
  end

  methods (Access=protected)
    
    function setDefaultValues(app)
      app.VariableName.Value = '';
      app.ShapeDropdown.Value = '';
      app.ViscositySpinner.Value = 1;
      app.updateCheckButton.Value = true;
      app.UITable.Data = repmat({''}, 6, 6);
    end
    
    function code = generateCode(app)
      code = {};  % TODO
    end
    
    function createMainComponents(app)
      % Create application specific components
      
      % Disable window resize
      app.UIFigure.Resize = 'off';
      
      % Create grid
      app.MainGrid = uigridlayout(app.UIFigure);
      app.MainGrid.ColumnWidth = {200, '1x'};
      app.MainGrid.RowHeight = {32, 32, 32, '1x', 32};
      app.MainGrid.RowSpacing = 0;
      app.MainGrid.ColumnSpacing = 20;
      
      % Variable name
      app.VariableName = ott.ui.support.OutputVariableEntry(app.MainGrid);
      app.VariableName.Layout.Row = 1;
      app.VariableName.Layout.Column = 1;
      
      % Shape selector
      app.ShapeDropdown = ott.ui.support.VariableDropdown(app.MainGrid, ...
          'label', 'Shape', 'filter', 'ott.shape.Shape');
      app.ShapeDropdown.Layout.Row = 2;
      app.ShapeDropdown.Layout.Column = 1;
      app.registerRefreshInput(app.ShapeDropdown);
        
      % Viscosity
      app.ViscositySpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Viscosity');
      app.ViscositySpinner.Layout.Row = 3;
      app.ViscositySpinner.Layout.Column = 1;
      
      % Auto-update widget
      app.updateCheckButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.updateCheckButton.Layout.Row = 5;
      app.updateCheckButton.Layout.Column = 1;
      
      % Create UITable
      cwidth = 40;
      app.UITable = uitable(app.MainGrid);
      app.UITable.Layout.Row = [1, 5];
      app.UITable.Layout.Column = 2;
      app.UITable.ColumnName = {'x', 'y', 'z', 'Rx', 'Ry', 'Rz'};
      app.UITable.RowName = {'x', 'y', 'z', 'Rx', 'Ry', 'Rz'};
      app.UITable.ColumnWidth = repmat({cwidth}, 1, 6);
      
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