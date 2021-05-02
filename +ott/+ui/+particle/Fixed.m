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

    aboutText = ['Generate a particle representation.' ...
      ' Combines the optical, fluidic, shape and other particle ' ...
      ' related properties into a single entity.'];
    
    helpText = {ott.ui.particle.Fixed.aboutText, ...
      '', ...
      ['Internally this applicatio creates a new ott.particle.Fixed' ...
      ' entitiy with parameters:'], ...
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
    ShapeDropdown       ott.ui.support.VariableDropdown
    DragDropdown        ott.ui.support.VariableDropdown
    TmatrixDropdown     ott.ui.support.VariableDropdown
    TinternalDropdown   ott.ui.support.VariableDropdown
    MassSpinner         ott.ui.support.LabeledSpinner
  end
  
  methods (Access=protected)
    function particle = generateParticle(app)
      % Generate the particle
      particle = ott.particle.Fixed(app.ShapeDropdown.Value, ...
          app.DragDropdown.Value, app.TmatrixDropdown.Value, ...
          'tinternal', app.TinternalDropdown.Value, ...
          'mass', app.MassSpinner.Value);
    end
    
    function code = generateCode(app)
      % Generate code for the user
      code = {};
      code{end+1} = 'particle = ott.particle.Fixed(app.ShapeDropdown.Value, ...';
      code{end+1} = '    app.DragDropdown.Value, app.TmatrixDropdown.Value, ...';
      code{end+1} = '    ''tinternal'', app.TinternalDropdown.Value, ...';
      code{end+1} = '    ''mass'', app.MassSpinner.Value);';
    end
    
    function setDefaultValues(app)
      % Set App-specific default values
      app.ShapeDropdown.Value = '';
      app.DragDropdown.Value = '';
      app.TmatrixDropdown.Value = '';
      app.TinternalDropdown.Value = '';
      app.MassSpinner.Value = 1e-6;
      
      setDefaultValues@ott.ui.particle.NewParticleBase(app);
    end
    
    function createMainComponents(app)
      % Create app specific components
      
      % Call base
      createMainComponents@ott.ui.particle.NewParticleBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = {32, 32, 32};
      
      % Shape
      app.ShapeDropdown = ott.ui.support.VariableDropdown(app.ExtraGrid, ...
          'label', 'Shape');
      app.ShapeDropdown.Layout.Row = 1;
      app.ShapeDropdown.Layout.Column = 1;
      app.registerRefreshInput(app.ShapeDropdown);
        
      % Drag
      app.DragDropdown = ott.ui.support.VariableDropdown(app.ExtraGrid, ...
          'label', 'Drag');
      app.DragDropdown.Layout.Row = 2;
      app.DragDropdown.Layout.Column = 1;
      app.registerRefreshInput(app.DragDropdown);
      
      % Mass
      app.MassSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Mass');
      app.MassSpinner.Layout.Row = 3;
      app.MassSpinner.Layout.Column = 1;
      
      % T-matrix
      app.TmatrixDropdown = ott.ui.support.VariableDropdown(app.ExtraGrid, ...
          'label', 'T-matrix');
      app.TmatrixDropdown.Layout.Row = 1;
      app.TmatrixDropdown.Layout.Column = 2;
      app.registerRefreshInput(app.TmatrixDropdown);
      
      % T-internal
      app.TinternalDropdown = ott.ui.support.VariableDropdown(app.ExtraGrid, ...
          'label', 'Internal T-matrix');
      app.TinternalDropdown.Layout.Row = 2;
      app.TinternalDropdown.Layout.Column = 2;
      app.registerRefreshInput(app.TinternalDropdown);
      
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
