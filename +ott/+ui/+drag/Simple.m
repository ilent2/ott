classdef Simple < ott.ui.support.AppTopLevel ...
    & ott.ui.support.RefreshInputsMenu ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.AppProducer
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
    ShapeDropDown         ott.ui.support.VariableDropDown
    ViscositySpinner      ott.ui.support.LabeledSpinner
    
    % Right column
    UITable                             matlab.ui.control.Table
  end

  methods (Access=protected)
    
    function setDefaultValues(app)
      app.VariableName.Value = '';
      app.ShapeDropDown.Value = '';
      app.ViscositySpinner.Value = 1.7e-3;
      app.UpdateButton.AutoUpdate = true;
      app.UITable.Data = repmat({''}, 6, 6);
    end
    
    function code = generateCode(app)
      % Generate code
      code = {};
      code{end+1} = ['viscosity = ' num2str(app.ViscositySpinner.Value) ';'];
      code{end+1} = 'drag = ott.drag.Stokes.FromShape(...';
      code{end+1} = '  shape, viscosity);';
    end
    
    function data = generateData(app)
      % Generate data
      data = ott.drag.Stokes.FromShape(...
        app.ShapeDropDown.Variable, ...
        app.ViscositySpinner.Value);
    end
    
    function updatePreview(app)
      % Update the preview table
      
      % Get data
      data = app.Data.forward;
      
      % Get string represeentation
      data = arrayfun(@(x)sprintf('%e', x), data, 'UniformOutput', false);
      
      % Display data
      app.UITable.Data = data;
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
      app.ShapeDropDown = ott.ui.support.VariableDropDown(app.MainGrid, ...
          'label', 'Shape', 'filter', 'ott.shape.Shape');
      app.ShapeDropDown.Layout.Row = 2;
      app.ShapeDropDown.Layout.Column = 1;
      app.ShapeDropDown.ValueChangedFcn = @(~,~) app.updateParametersCb();
      app.registerRefreshInput(app.ShapeDropDown);
        
      % Viscosity
      app.ViscositySpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
          'label', 'Viscosity');
      app.ViscositySpinner.Layout.Row = 3;
      app.ViscositySpinner.Layout.Column = 1;
      app.ViscositySpinner.Step = 1e-4;
      app.ViscositySpinner.Limits = [0, Inf];
      app.ViscositySpinner.LowerLimitInclusive = 'off';
      app.ViscositySpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Auto-update widget
      app.UpdateButton = ott.ui.support.UpdateCheckButton(app.MainGrid);
      app.UpdateButton.Layout.Row = 5;
      app.UpdateButton.Layout.Column = 1;
      
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
      
      % Construct widgets first
      app = app@ott.ui.support.AppTopLevel();
      
      % Then connect producer callbacks
      app = app@ott.ui.support.AppProducer();
      
      if nargout == 0
        clear app;
      end
    end
  end
end