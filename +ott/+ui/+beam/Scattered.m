classdef Scattered < ott.ui.beam.NewBeamBase ...
    & ott.ui.support.RefreshInputsMenu
% GUI to generate a scattered beam from a particle and beam input.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Scattered';

    nameText = 'Scatter Beam';

    aboutText = ['Generate a scattered beam from a particle and incident' ...
      'beam instance.'];
    
    helpText = {ott.ui.beam.Scattered.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.Scattered.nameText;
  end
  
  properties (Access=public)
    IncidentBeamDropDown      ott.ui.support.VariableDropDown
    ParticleDropDown          ott.ui.support.VariableDropDown
    ParticlePanel             matlab.ui.container.Panel
    ParticlePosition          ott.ui.support.XyzSpinners
    ParticleRotation          ott.ui.support.XyzSpinners
  end
  
  methods (Access=protected)
    
    function setDefaultValues(app)
      
      % Add polarisation default value
      app.IncidentBeamDropDown.Value = '';
      app.ParticleDropDown.Value = '';
      app.ParticlePosition.Value = [0,0,0];
      app.ParticleRotation.Value = [0,0,0];
      
      % Call remaining defaults
      setDefaultValues@ott.ui.beam.NewBeamBase(app);
    end
    
    function code = generateCode(app)
      code = {};   % TODO
    end
    
    function beam = generateData(app)
      try
        beam = app.IncidentBeamDropDown.Value.scatter(...
            app.ParticleDropDown.Value, ...
            'position', app.ParticlePosition.Value, ...
            'rotation', app.ParticleRotation.Value);
        beam.position = app.PositionXyzSpinner.Value;
        beam.rotation = app.RotationXyzSpinner.Value;
      catch
        beam = [];
      end
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = {32, 32, 90, '1x'};
      
      % Hide obselete rows of main grid
      app.MainGrid.RowHeight{app.WavelengthSpinner.Layout.Row} = 0;
      app.MainGrid.RowHeight{app.IndexSpinner.Layout.Row} = 0;

      % Incident beam drop down
      app.IncidentBeamDropDown = ott.ui.support.VariableDropDown(...
        app.ExtraGrid, 'label', 'Incident Beam');
      app.IncidentBeamDropDown.Layout.Row = 1;
      app.registerRefreshInput(app.IncidentBeamDropDown);
      
      % Particle drop down
      app.ParticleDropDown = ott.ui.support.VariableDropDown(...
        app.ExtraGrid, 'label', 'Particle');
      app.ParticleDropDown.Layout.Row = 2;
      app.registerRefreshInput(app.ParticleDropDown);
      
      % Particle properties panel
      app.ParticlePanel = uipanel(app.ExtraGrid);
      app.ParticlePanel.Title = 'Particle Properties';
      app.ParticlePanel.Layout.Row = 3;
      
      % Grid
      grid = uigridlayout(app.ParticlePanel);
      grid.Padding = [5, 5, 5, 5];
      grid.ColumnWidth = {'1x'};
      grid.ColumnSpacing = 1;
      grid.RowSpacing = 1;
      grid.RowHeight = {32, 32};
      
      % Particle position entry
      app.ParticlePosition = ott.ui.support.XyzSpinners(...
          grid, 'label', 'Pos.');
      app.ParticlePosition.Layout.Row = 1;
      app.ParticlePosition.Layout.Column = 1;
      app.ParticlePosition.Step = 1e-7;
      app.ParticlePosition.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Particle rotation entry
      app.ParticleRotation = ott.ui.support.XyzSpinners(...
          grid, 'label', 'Rot.');
      app.ParticleRotation.Layout.Row = 2;
      app.ParticleRotation.Layout.Column = 1;
      app.ParticleRotation.Step = 10;
      app.ParticleRotation.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
    end
  end
  
  methods (Access=public)
    function app=Scattered()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.NewBeamBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end