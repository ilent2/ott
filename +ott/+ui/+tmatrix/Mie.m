classdef Mie < ott.ui.tmatrix.NewTmatrixBase
% Generate a T-matrix cotaining the Mie coefficients.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Mie';

    nameText = 'Mie T-matrix';

    aboutText = ['Generate T-matrix for a spherical particle using ' ...
      'the Mie coefficients.'];
    
    helpText = {ott.ui.tmatrix.Mie.aboutText, ...
      ''};
    
    windowName = ott.ui.tmatrix.Mie.nameText;
  end
  
  properties (Access=public)
    InputDropDown         ott.ui.support.LabeledDropDown
    RadiusSpinner         ott.ui.support.LabeledSpinner
  end
  
  methods (Access=protected)
    function  setDefaultValues(app)
      
      % Set this widgets defaults
      app.InputDropDown.Value = 'Shape Max Radius';
      app.RadiusSpinner.Value = 1e-6;
      
      % Change widget status based on InputDropDown
      app.updateVisibleWidgets();
      
      % Defer to base
      setDefaultValues@ott.ui.tmatrix.NewTmatrixBase(app);
    end
    
    function data = generateData(app)
      % Generate shape
      
      switch app.InputDropDown.Value
        case 'Radius'
          radius = app.RadiusSpinner.Value;
        case 'Shape Max Radius'
          radius = app.ShapeName.Variable.maxRadius;
        case 'Shape Volume'
          radius = ((3/4/pi) * app.ShapeName.Variable.volume).^(1/3);
        otherwise
          error('Internal error');
      end
      
      data = ott.tmatrix.Mie(...
        'radius', radius ./ app.WavelengthSpinner.Value, ...
        'index_relative', app.RelativeIndexSpinner.Value);
    end
    
    function code = generateCode(app)
      code = {};
      
      switch app.InputDropDown.Value
        case 'Radius'
          code{end+1} = ['radius = ' num2str(app.RadiusSpinner.Value) ';'];
        case 'Shape Max Radius'
          code{end+1} = ['radius = ' app.ShapeName.Value '.maxRadius;'];
        case 'Shape Volume'
          code{end+1} = ['radius = ((3/4/pi) * ' app.ShapeName.Value '.volume).^(1/3);'];
        otherwise
          error('Internal error');
      end
      
      code{end+1} = ['wavelength = ' num2str(app.WavelengthSpinner.Value) ';'];
      code{end+1} = ['index_relative = ' num2str(app.RelativeIndexSpinner.Value) ';'];
      code{end+1} = '';
      code{end+1} = 'tmatrix = ott.tmatrix.Mie(...';
      code{end+1} = '  ''radius'', radius ./ wavelength, ...';
      code{end+1} = '  ''index_relative'', index_relative);';
    end
    
    function updateVisibleWidgets(app)
      % Update widget enable status
      if strcmpi(app.InputDropDown.Value, 'Radius')
        app.RadiusSpinner.Enable = true;
        app.ShapeName.Enable = false;
        app.MainGrid.RowHeight{3} = app.MainGrid.RowHeight{1};
        app.MainGrid.RowHeight{4} = 0;
      else
        app.RadiusSpinner.Enable = false;
        app.ShapeName.Enable = true;
        app.MainGrid.RowHeight{3} = 0;
        app.MainGrid.RowHeight{4} = app.MainGrid.RowHeight{1};
      end
    end
    
    function inputTypeChangedCb(app)
      
      % Update visible widgets
      app.updateVisibleWidgets();
      
      % Run auto-update
      app.updateParametersCb();
    end
    
    function createMainComponents(app)
      % Create components for gui
      
      % Create base components
      createMainComponents@ott.ui.tmatrix.NewTmatrixBase(app);
      
      % Increase the main grid size and shift everything down by one
      % So we can put the Shape Type input after the Output field
      ott.ui.support.insertGridRow(app.MainGrid, app.MainGrid.RowHeight{1}, 2, 2);
      
      % Add Input type selector bellow output selector
      app.InputDropDown = ott.ui.support.LabeledDropDown(app.MainGrid, ...
        'label', 'Input Type', 'wwidth', app.wwidth);
      app.InputDropDown.Layout.Row = 2;
      app.InputDropDown.Layout.Column = 1;
      app.InputDropDown.Items = {'Shape Max Radius', 'Shape Volume', 'Radius'};
      app.InputDropDown.ValueChangedFcn = @(~,~) app.inputTypeChangedCb();
      
      % Radius entry
      app.RadiusSpinner = ott.ui.support.LabeledSpinner(app.MainGrid, ...
        'label', 'Radius', 'wwidth', app.wwidth);
      app.RadiusSpinner.Layout.Row = 3;
      app.RadiusSpinner.Layout.Column = 1;
      app.RadiusSpinner.Step = 1e-7;
      app.RadiusSpinner.Limits = [0, Inf];
      app.RadiusSpinner.LowerLimitInclusive = 'off';
      app.RadiusSpinner.ValueChangedFcn = @(~,~) app.inputTypeChangedCb();
      
    end
  end
  
  methods (Access=public)
    function app=Mie()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.NewTmatrixBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end