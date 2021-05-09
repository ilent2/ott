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
  end
  
  methods (Access=protected)
    
    function setDefaultValues(app)
      
      % Add polarisation default value
      app.IncidentBeamDropDown.Value = '';
      app.ParticleDropDown.Value = '';
      
      % Call remaining defaults
      setDefaultValues@ott.ui.beam.NewBeamBase(app);
    end
    
    function code = generateCode(app)
      code = {};
    end
    
    function beam = generateData(app)
      try
        beam = app.IncidentBeamDropDown.Value.scatter(...
            app.ParticleDropDown.Value, ...
            'position', app.TranslationXyzSpinner.Value, ...
            'rotation', app.RotationXyzSpinner.Value);
      catch
        beam = [];
      end
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);
      
      % Hide obselete rows of main grid
      app.MainGrid.RowHeight{app.WavelengthSpinner.Layout.Row} = 0;
      app.MainGrid.RowHeight{app.IndexSpinner.Layout.Row} = 0;
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = {32, 32, '1x'};

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