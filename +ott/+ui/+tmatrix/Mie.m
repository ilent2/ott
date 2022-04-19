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
        case 'Sphere'
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
        'index_relative', app.RelativeIndexSpinner);
    end
    
    function code = generateCode(app)
      code = {}; % TODO
      
      switch app.InputDropDown.Value
        case 'Sphere'
        case 'Shape Max Radius'
        case 'Shape Volume'
        otherwise
          error('Internal error');
      end
    end
    
    function updateVisibleWidgets(app)
      % Update widget enable status
      if strcmpi(app.InputDropDown.Value, 'Sphere')
        app.RadiusSpinner.Enable = true;
        app.ShapeName.Enable = false;
      else
        app.RadiusSpinner.Enable = false;
        app.ShapeName.Enable = true;
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
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 3);
      app.ExtraGrid.RowHeight(end) = {'1x'};
      
      % Input type selector
      app.InputDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
        'label', 'Input', 'wwidth', app.wwidth);
      app.InputDropDown.Layout.Row = 1;
      app.InputDropDown.Layout.Column = 1;
      app.InputDropDown.Items = {'Shape Max Radius', 'Shape Volume', 'Sphere'};
      app.InputDropDown.ValueChangedFcn = @(~,~) app.inputTypeChangedCb();
      
      % Radius entry
      app.RadiusSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Radius', 'wwidth', app.wwidth);
      app.RadiusSpinner.Layout.Row = 2;
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