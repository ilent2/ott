classdef FromShape < ott.ui.particle.NewParticleBase ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Create a simple particle using the particle.Fixed.FromShape method.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Progress bar (too hard for now, add later)

  properties (Constant)
    cnameText = 'FromShape';

    nameText = 'Generate Particle from Shape';

    aboutText = ['Generate a simple particle from a shape. ' ...
      'Calculates the optical, fluidic, shape and other particle ' ...
      'related properties into a single entity.'];
    
    helpText = {ott.ui.particle.FromShape.aboutText, ...
      '', ...
      ['Internally this application uses the ott.particle.Fixed.FromShape' ...
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
      'Mass (kg) -- Mass of the particle.'};
    
    windowName = ott.ui.particle.FromShape.nameText;
    windowSize = [500, 260];
  end
  
  properties (Access=public)
    ShapeDropDown       ott.ui.support.VariableDropDown
    ParticleIrSpinner   ott.ui.support.LabeledSpinner
    MediumIrSpinner     ott.ui.support.LabeledSpinner
    WavelengthSpinner   ott.ui.support.LabeledSpinner
    ViscositySpinner    ott.ui.support.LabeledSpinner
    InternalCheckbox    matlab.ui.control.CheckBox
    MassSpinner         ott.ui.support.LabeledSpinner
  end
  
  methods (Access=protected)
    function data = generateData(app)
      % Generate particle
      
      data = ott.particle.Fixed.FromShape(app.ShapeDropDown.Variable, ...
        'index_particle', app.ParticleIrSpinner.Value, ...
        'index_medium', app.MediumIrSpinner.Value, ...
        'wavelength0', app.WavelengthSpinner.Value, ...
        'viscosity', app.ViscositySpinner.Value, ...
        'internal', app.InternalCheckbox.Value, ...
        'mass', app.MassSpinner.Value);
    end
    
    function code = generateCode(app)
      % Generate code
      
      % Get default shape name
      shape_name = app.ShapeDropDown.Value;
      if isempty(app.ShapeDropDown.Value)
        shape_name = 'shape';
      end
      
      code = {};
      code{end+1} = ['particle = ott.particle.Fixed.FromShape(', shape_name ', ...'];
      code{end+1} = ['  ''index_particle'', ', num2str(app.ParticleIrSpinner.Value), ', ...'];
      code{end+1} = ['  ''index_medium'', ', num2str(app.MediumIrSpinner.Value), ', ...'];
      code{end+1} = ['  ''wavelength0'', ', num2str(app.WavelengthSpinner.Value), ', ...'];
      code{end+1} = ['  ''viscosity'', ', num2str(app.ViscositySpinner.Value), ', ...'];
      if app.InternalCheckbox.Value
        code{end+1} = '  ''internal'', true, ...';
      else
        code{end+1} = '  ''internal'', false, ...';
      end
      code{end+1} = ['  ''mass'', ', num2str(app.MassSpinner.Value), ');'];
    end
    
    function setDefaultValues(app)
      % Set App-specific default values
      
      app.ShapeDropDown.Value = '';
      app.ParticleIrSpinner.Value = 1.2;
      app.MediumIrSpinner.Value = 1.0;
      app.WavelengthSpinner.Value = 1e-6;
      app.ViscositySpinner.Value = 8.9e-4;
      app.MassSpinner.Value = 1e-6;
      app.InternalCheckbox.Value = false;
      
      setDefaultValues@ott.ui.particle.NewParticleBase(app);
    end
    
    function createMainComponents(app)
      % Create app specific components
      
      % Call base
      createMainComponents@ott.ui.particle.NewParticleBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = repmat({30}, 1, 4);
      
      % Hide progress bar
      % TODO: This GUI should have a progress bar, but it should work for
      % the particle T-matrix calculation, drag calculation and anything
      % else.  Not sure how do to this at the moment.
      app.UpdateButton.Gauge.Visible = false;
      app.UpdateButton.Button.Layout.Column = 2;
      
      % Shape
      app.ShapeDropDown = ott.ui.support.VariableDropDown(app.ExtraGrid, ...
          'label', 'Shape', 'filter', 'ott.shape.Shape');
      app.ShapeDropDown.Layout.Row = 1;
      app.ShapeDropDown.Layout.Column = 1;
      app.registerRefreshInput(app.ShapeDropDown);
      
      % Particle IR spinner
      app.ParticleIrSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Particle Index');
      app.ParticleIrSpinner.Step = 0.1;
      app.ParticleIrSpinner.Limits = [0, Inf];
      app.ParticleIrSpinner.LowerLimitInclusive = 'off';
      app.ParticleIrSpinner.Layout.Row = 2;
      app.ParticleIrSpinner.Layout.Column = 1;
        
      % Medium IR spinner
      app.MediumIrSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Medium Index');
      app.MediumIrSpinner.Step = 0.1;
      app.MediumIrSpinner.Limits = [0, Inf];
      app.MediumIrSpinner.LowerLimitInclusive = 'off';
      app.MediumIrSpinner.Layout.Row = 3;
      app.MediumIrSpinner.Layout.Column = 1;
      
      % Internal checkbox
      app.InternalCheckbox = uicheckbox(app.ExtraGrid);
      app.InternalCheckbox.Text = 'Internal';
      app.InternalCheckbox.Layout.Row = 4;
      app.InternalCheckbox.Layout.Column = 1;
        
      % Wavelength
      app.WavelengthSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Wavelength');
      app.WavelengthSpinner.Step = 1e-7;
      app.WavelengthSpinner.Limits = [0, Inf];
      app.WavelengthSpinner.LowerLimitInclusive = 'off';
      app.WavelengthSpinner.Layout.Row = 1;
      app.WavelengthSpinner.Layout.Column = 2;
        
      % Viscosity
      app.ViscositySpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Viscosity');
      app.ViscositySpinner.Step = 1e-5;
      app.ViscositySpinner.Limits = [0, Inf];
      app.ViscositySpinner.LowerLimitInclusive = 'off';
      app.ViscositySpinner.Layout.Row = 2;
      app.ViscositySpinner.Layout.Column = 2;
        
      % Mass
      app.MassSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Mass');
      app.MassSpinner.Step = 1e-7;
      app.MassSpinner.Limits = [0, Inf];
      app.MassSpinner.LowerLimitInclusive = 'off';
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