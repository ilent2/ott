classdef Fixed < ott.ui.particle.NewParticleBase ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Create a particle with fixed properties.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Fixed';

    nameText = 'Generate Simple Particle';

    aboutText = ['Generate a particle representation. ' ...
      'Combines the optical, fluidic, shape and other particle ' ...
      'related properties into a single entity.'];
    
    helpText = {ott.ui.particle.Fixed.aboutText, ...
      '', ...
      ['Internally this application creates a new ott.particle.Fixed' ...
      ' entity with parameters:'], ...
      '', ...
      'Shape -- Geometric shape.', ...
      '', ...
      'Drag -- Particle drag description.', ...
      '', ...
      'T-matrix -- External T-matrix.', ...
      '', ...
      'T-internal -- Internal T-matrix.', ...
      '', ...
      'Mass -- Mass of the particle.'};
    
    windowName = ott.ui.particle.Fixed.nameText;
    windowSize = [500, 230];
  end
  
  properties
    ShapeDropDown       ott.ui.support.VariableDropDown
    DragDropDown        ott.ui.support.VariableDropDown
    TmatrixDropDown     ott.ui.support.VariableDropDown
    TinternalDropDown   ott.ui.support.VariableDropDown
    MassSpinner         ott.ui.support.LabeledSpinner
  end
  
  methods (Access=protected)
    function particle = generateData(app)
      % Generate the particle
      particle = ott.particle.Fixed(app.ShapeDropDown.Variable, ...
          app.DragDropDown.Variable, app.TmatrixDropDown.Variable, ...
          'tinternal', app.TinternalDropDown.Variable, ...
          'mass', app.MassSpinner.Value);
    end
    
    function code = generateCode(app)
      % Generate code for the user
      code = {};
      code{end+1} = 'particle = ott.particle.Fixed(app.ShapeDropDown.Value, ...';
      code{end+1} = '    app.DragDropDown.Value, app.TmatrixDropDown.Value, ...';
      code{end+1} = '    ''tinternal'', app.TinternalDropDown.Value, ...';
      code{end+1} = '    ''mass'', app.MassSpinner.Value);';
    end
    
    function setDefaultValues(app)
      % Set App-specific default values
      app.ShapeDropDown.Value = '';
      app.DragDropDown.Value = '';
      app.TmatrixDropDown.Value = '';
      app.TinternalDropDown.Value = '';
      app.MassSpinner.Value = 1e-6;
      
      setDefaultValues@ott.ui.particle.NewParticleBase(app);
    end
    
    function createMainComponents(app)
      % Create app specific components
      
      % Call base
      createMainComponents@ott.ui.particle.NewParticleBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = {32, 32, 32};
      
      % Hide progress bar
      app.UpdateButton.Gauge.Visible = false;
      app.UpdateButton.Button.Layout.Column = 2;
      
      % Shape
      app.ShapeDropDown = ott.ui.support.VariableDropDown(app.ExtraGrid, ...
          'label', 'Shape', 'filter', 'ott.shape.Shape');
      app.ShapeDropDown.Layout.Row = 1;
      app.ShapeDropDown.Layout.Column = 1;
      app.registerRefreshInput(app.ShapeDropDown);
        
      % Drag
      app.DragDropDown = ott.ui.support.VariableDropDown(app.ExtraGrid, ...
          'label', 'Drag', 'filter', 'ott.drag.Stokes');
      app.DragDropDown.Layout.Row = 2;
      app.DragDropDown.Layout.Column = 1;
      app.registerRefreshInput(app.DragDropDown);
      
      % Mass
      app.MassSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Mass');
      app.MassSpinner.Layout.Row = 3;
      app.MassSpinner.Layout.Column = 1;
      app.MassSpinner.Step = 0.5e-6;
      
      % T-matrix
      app.TmatrixDropDown = ott.ui.support.VariableDropDown(app.ExtraGrid, ...
          'label', 'T-matrix', 'filter', 'ott.tmatrix.Tmatrix');
      app.TmatrixDropDown.Layout.Row = 1;
      app.TmatrixDropDown.Layout.Column = 2;
      app.registerRefreshInput(app.TmatrixDropDown);
      
      % T-internal
      app.TinternalDropDown = ott.ui.support.VariableDropDown(app.ExtraGrid, ...
          'label', 'Internal T-matrix', 'filter', 'ott.tmatrix.Tmatrix');
      app.TinternalDropDown.Layout.Row = 2;
      app.TinternalDropDown.Layout.Column = 2;
      app.registerRefreshInput(app.TinternalDropDown);
      
    end
  end
  
  methods (Access=public)
    function app=Fixed()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.particle.NewParticleBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
end
