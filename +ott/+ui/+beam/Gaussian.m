classdef Gaussian < ott.ui.beam.NewBeamBase ...
    & ott.ui.support.GenerateCodeMenu
% Generate a simple beam representation and visualise.
%
% Supported beams:
%   - Gaussian
%   - Laguerre-Gaussian
%   - Hermite-Gaussian
%
% Some of these beams might move to their own interface in a future
% release.
%
% This GUI can be launched from the launcher Beam -> Simple or by running:
%
%   ott.ui.beam.Simple

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Do we want a NA warning?

  properties (Constant)
    cnameText = 'Gaussian';

    nameText = 'Create a Gaussian Beam';

    aboutText = ['Generate a Gaussian or related beam.  Can be used to' ...
      ' create Gaussian, Laguerre-Gaussian, Hermite-Gaussian.'];
    
    helpText = {ott.ui.beam.Gaussian.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.Gaussian.nameText;
  end
  
  properties (Access=public)
    BeamTypeDropDown            ott.ui.support.LabeledDropDown
    LmodeSpinner                ott.ui.support.LabeledSpinner
    PmodeSpinner                ott.ui.support.LabeledSpinner
    MmodeSpinner                ott.ui.support.LabeledSpinner
    NmodeSpinner                ott.ui.support.LabeledSpinner
    PowerSpinner                ott.ui.support.LabeledSpinner
    NaSpinner                   ott.ui.support.LabeledSpinner
    PolarisationEntry           ott.ui.support.JonesPolarisationEntry
  end
  
  methods (Access=protected)
    
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function data = generateData(app)
      data = []; % TODO
    end
    
    function TypeDropDownValueChanged(app, evt)
      
      % Update visible shape widgets
      app.UpdateVisibleTypeWidgets();
      
      % Continue processing changed value
      app.updateParametersCb(evt);
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);
      
      % Create grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 15);
      app.ExtraGrid.RowHeight{end} = '1x';
      app.ExtraGrid.RowSpacing = 1;
        
      % Beam type dropdown
      app.BeamTypeDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid);
      app.BeamTypeDropDown.Layout.Row = 2;
      app.BeamTypeDropDown.Layout.Column = 1;
      app.BeamTypeDropDown.ValueChangedFcn = createCallbackFcn(app, ...
          @typeChangedCb, true);
        
      % Lmode spinner
      app.LmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'L mode');
      app.LmodeSpinner.Layout.Row = 3;
      app.LmodeSpinner.Layout.Column = 1;
      app.LmodeSpinner.Step = 1;
      app.LmodeSpinner.Limits = [-Inf,Inf];
      app.LmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Pmode spinner
      app.PmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'P mode');
      app.PmodeSpinner.Layout.Row = 4;
      app.PmodeSpinner.Layout.Column = 1;
      app.PmodeSpinner.Step = 1;
      app.PmodeSpinner.Limits = [0,Inf];
      app.PmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Mmode spinner
      app.MmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'M mode');
      app.MmodeSpinner.Layout.Row = 5;
      app.MmodeSpinner.Layout.Column = 1;
      app.MmodeSpinner.Step = 1;
      app.MmodeSpinner.Limits = [0,Inf];
      app.MmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Nmode spinner
      app.NmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'N mode');
      app.NmodeSpinner.Layout.Row = 6;
      app.NmodeSpinner.Layout.Column = 1;
      app.NmodeSpinner.Step = 1;
      app.NmodeSpinner.Limits = [0,Inf];
      app.NmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
        
      % NA spinner
      app.NaSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Numerical Aperture');
      app.NaSpinner.Layout.Row = 9;
      app.NaSpinner.Layout.Column = 1;
      app.NaSpinner.Step = 0.1;
      app.NaSpinner.LowerLimitInclusive = 'off';
      app.NaSpinner.Limits = [0.0, Inf];
      app.NaSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
        
      % Polarisation jones vector entry
      app.PolarisationEntry = ott.ui.support.JonesPolarisationEntry(...
          app.ExtraGrid);
      app.PolarisationEntry.Layout.Row = 10;
      app.PolarisationEntry.Layout.Column = 1;
      app.PolarisationEntry.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Power spinner  
      app.PowerSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Power');
      app.PowerSpinner.Layout.Row = 11;
      app.PowerSpinner.Layout.Column = 1;
      app.PowerSpinner.Step = 0.1;
      app.PowerSpinner.LowerLimitInclusive = 'off';
      app.PowerSpinner.Limits = [0.0, Inf];
      app.PowerSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();

    end
  end
  
  methods (Access=public)
    function app=Gaussian()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.NewBeamBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end