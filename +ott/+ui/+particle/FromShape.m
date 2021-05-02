classdef FromShape < ott.ui.particle.NewParticleBase ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Create a simple particle using the particle.Fixed.FromShape method.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'FromShape';

    nameText = 'Generate Particle from Shape';

    aboutText = ['Generate a simple particle from a shape.' ...
      ' Calculates the optical, fluidic, shape and other particle ' ...
      ' related properties into a single entity.'];
    
    helpText = {ott.ui.particle.FromShape.aboutText, ...
      '', ...
      ['Internally this applicatio uses the ott.particle.Fixed.FromShape' ...
      ' method to create the particle.  Options are:'], ...
      '', ...
      'Shape -- The geometric shape describing the object geometry.', ...
      '', ...
      'Particle Index -- Refractive index of the particle.', ...
      '', ...
      'Medium Index -- Refractive index of the medium.', ...
      '', ...
      'Wavelength Medium -- Wavelength in the medium.', ...
      '', ...
      'Viscosity (Ns/m2) -- Viscosity of the medium.', ...
      '', ...
      'Calculate Internal -- If true, also calculates internal T-matrix.', ...
      '', ...
      'Mass -- Mass of the particle.'};
    
    windowName = ott.ui.particle.FromShape.nameText;
    windowSize = [500, 230];
  end
  
  properties (Access=public)
    ShapeDropdown       ott.ui.support.VariableDropdown
    ParticleIrSpinner    ott.ui.support.LabeledSpinner
    MediumIrSpinner    ott.ui.support.LabeledSpinner
    WavelengthSpinner    ott.ui.support.LabeledSpinner
    ViscositySpinner    ott.ui.support.LabeledSpinner
    InternalCheckbox    matlab.ui.control.CheckBox
    MassSpinner         ott.ui.support.LabeledSpinner
  end
  
  methods (Access=protected)
    function particle = generateParticle(app)
      % TODO
      particle = [];
    end
    
    function code = generateCode(app)
      % TODO
      code = {};
    end
    
    function setDefaultValues(app)
      % Set App-specific default values
      
      app.ShapeDropdown.Value = '';
      app.ParticleIrSpinner.Value = 1.2;
      app.MediumIrSpinner.Value = 1.0;
      app.WavelengthSpinner.Value = 1e-6;
      app.ViscositySpinner.Value = 1e-3;
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
      
      % Particle IR spinner
      app.ParticleIrSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Particle Index');
      app.ParticleIrSpinner.Layout.Row = 2;
      app.ParticleIrSpinner.Layout.Column = 1;
        
      % Medium IR spinner
      app.MediumIrSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Medium Index');
      app.MediumIrSpinner.Layout.Row = 3;
      app.MediumIrSpinner.Layout.Column = 1;
        
      % Wavelength
      app.WavelengthSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Wavelength');
      app.WavelengthSpinner.Layout.Row = 1;
      app.WavelengthSpinner.Layout.Column = 2;
        
      % Viscosity
      app.ViscositySpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Viscosity');
      app.ViscositySpinner.Layout.Row = 2;
      app.ViscositySpinner.Layout.Column = 2;
        
      % Mass
      app.MassSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Mass');
      app.MassSpinner.Layout.Row = 3;
      app.MassSpinner.Layout.Column = 2;
      
    end
  end
  
  methods (Access=public)
    function app=FromShape()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.particle.NewParticleBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
end