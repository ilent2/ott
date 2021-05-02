classdef PlaneWave < ott.ui.beam.NewBeamBase
% Generate a plane wave beam.
%
% This GUI can be launched from the launcher Beam -> PlaneWave or with:
%
%   ott.ui.beam.PlaneWave

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'PlaneWave';

    nameText = 'Create Plane Wave Beam';

    aboutText = ['Generate a plane wave beam.'];
    
    helpText = {ott.ui.beam.PlaneWave.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PlaneWave.nameText;
  end
  
  properties (Access=public)
    PolarisationEntry     ott.ui.support.JonesPolarisationEntry
  end
  
  methods (Access=protected)
    
    function setDefaultValues(app)
      
      % Add polarisation default value
      app.PolarisationEntry.Value = [1, 1i];
      
      % Call remaining defaults
      setDefaultValues@ott.ui.beam.NewBeamBase(app);
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);
      
      % Hide obselete rows of main grid
      app.MainGrid.RowHeight{app.TranslationXyzSpinner.Layout.Row} = 0;
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = {32, '1x'};
      
      % Polarisation jones vector entry
      app.PolarisationEntry = ott.ui.support.JonesPolarisationEntry(...
          app.ExtraGrid);
      app.PolarisationEntry.Layout.Row = 1;
      app.PolarisationEntry.Layout.Column = 1;
      app.PolarisationEntry.ValueChangedFcn = createCallbackFcn(app, ...
          @valueChangedCb, true);
      
    end
    
    function beam = generateBeam(app)
      
      rot3 = app.RotationXyzSpinner.Value;
      rotation = ott.utils.rotx(rot3(1))*ott.utils.roty(rot3(2))*ott.utils.rotz(rot3(3));
      
      % Generate new beam
      beam = ott.beam.PlaneWave(...
          'polarisation', app.PolarisationEntry.Value, ...
          'index_medium', app.IndexSpinner.Value, ...
          'wavelength0', app.WavelengthSpinner.Value, ...
          'rotation', rotation);
    end
  end
  
  methods (Access=public)
    function app=PlaneWave()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.NewBeamBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
end
