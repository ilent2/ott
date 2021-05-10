classdef PmParaxial < ott.ui.beam.NewBeamBase ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Generate a beam using paraxial point matching

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'PmParaxial';

    nameText = 'Paraxial Pointmatched Beam';

    aboutText = ['Generates a beam using paraxial point matching from' ...
      ' a given E field profile.'];
    
    helpText = {ott.ui.beam.PmParaxial.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
  end
  
  properties (Access=protected)
    ProfileDropDown           ott.ui.support.VariableDropDown
    MappingDropDown           ott.ui.support.LabeledDropDown
    MaxAngleSpinner           ott.ui.support.LabeledSpinner
    PolarisationDropDown      ott.ui.support.LabeledDropDown
    PolarisationEntry         ott.ui.support.JonesPolarisationEntry
  end
  
  methods (Access=protected)
    
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function data = generateData(app)
      data = []; % TODO
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);
      
      % Configure grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 6);
      app.ExtraGrid.RowHeight(end) = {'1x'};
      
      % TODO: Hide extra components
      
      % Input field
      app.ProfileDropDown = ott.ui.support.VariableDropDown(app.ExtraGrid, ...
        'label', 'Field');
      app.registerRefreshInput(app.ProfileDropDown);
      
      % Mapping
      app.MappingDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
        'label', 'Mapping');
      app.MappingDropDown.Items = {'sin', 'tan', 'theta'};
      app.MappingDropDown.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Maximum angle
      app.MaxAngleSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Max. Angle');
      app.MaxAngleSpinner.Step = 10;
      app.MaxAngleSpinner.Limits = [0, 90];
      app.MaxAngleSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Polarisation drop down
      app.PolarisationDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
        'label', 'Polarisation');
      app.PolarisationDropDown.Items = {'Polar', 'Cartesian'};
      app.PolarisationDropDown.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Polarisation value
      app.PolarisationEntry = ott.ui.support.JonesPolarisationEntry(app.ExtraGrid, ...
        'label', 'Polarisation');
      app.PolarisationEntry.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
    end
  end
  
  methods (Access=public)
    function app=PmParaxial()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.NewBeamBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end