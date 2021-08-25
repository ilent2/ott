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
%   ott.ui.beam.Gaussian

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Do we want a NA warning?
% TODO: Auto-update checkbox doesn't work, always updates with DropDown.

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
      
      switch app.BeamTypeDropDown.Value
        case 'Gaussian'
        case 'Hermite'
        case 'Laguerre'
        otherwise
          error('Internal error');
      end
    end
    
    function data = generateData(app)
      data = [];
      
      wavelength = app.WavelengthSpinner.Value;
      index = app.IndexSpinner.Value;
      position = app.PositionXyzSpinner.Value;
      rotation = ott.utils.euler2rot(app.RotationXyzSpinner.Value.');
      NA = app.NaSpinner.Value;
      polarisation = app.PolarisationEntry.Value;
      power = app.PowerSpinner.Value;
      
      warning('wavlength not used yet'); % TODO
      
      switch app.BeamTypeDropDown.Value
        case 'Gaussian'
          data = ott.beam.Gaussian.FromNa(NA, 'power', power, ...
            'index_medium', index, 'polfield', polarisation, ...
            'position', position, 'rotation', rotation);
        case 'Hermite'
          data = ott.beam.HermiteGaussian.FromNa(NA, 'power', power, ...
            'mmode', app.MmodeSpinner.Value, 'nmode', app.NmodeSpinner.Value, ...
            'index_medium', index, 'polfield', polarisation, ...
            'position', position, 'rotation', rotation);
        case 'Laguerre'
          data = ott.beam.LaguerreGaussian.FromNa(NA, 'power', power, ...
            'lmode', app.LmodeSpinner.Value, 'pmode', app.PmodeSpinner.Value, ...
            'index_medium', index, 'polfield', polarisation, ...
            'position', position, 'rotation', rotation);
        otherwise
          error('Internal error');
      end
    end
    
    function setDefaultValues(app)
      
      app.BeamTypeDropDown.Value = 'Gaussian';
      app.LmodeSpinner.Value = 2;
      app.PmodeSpinner.Value = 0;
      app.MmodeSpinner.Value = 2;
      app.NmodeSpinner.Value = 2;
      app.PowerSpinner.Value = 2;
      app.NaSpinner.Value = 0.9;
      app.PolarisationEntry.Value = [1, 1i];
      
      % Show/hide widgets
      app.configureRowVisibility();
      
      % Do base work
      setDefaultValues@ott.ui.beam.NewBeamBase(app);
    end
    
    function configureRowVisibility(app)
      % Show/hide rows for different beam types
      
      rheight = 32;
      
      % Hide all widgets
      app.ExtraGrid.RowHeight(1:end) = {0};
      
      % Disable all widgets
      widgets = findobj(app.ExtraGrid.Children, '-property', 'Enable');
      set(widgets, 'Enable', 'off');
      
      % Make common widgets visible again
      app.ExtraGrid.RowHeight{app.BeamTypeDropDown.Layout.Row} = rheight;
      app.ExtraGrid.RowHeight{app.PowerSpinner.Layout.Row} = rheight;
      app.ExtraGrid.RowHeight{app.NaSpinner.Layout.Row} = rheight;
      app.ExtraGrid.RowHeight{app.PolarisationEntry.Layout.Row} = rheight;
      app.BeamTypeDropDown.Enable = 'on';
      app.NaSpinner.Enable = 'on';
      app.PolarisationEntry.Enable = 'on';
      app.PowerSpinner.Enable = 'on';
      
      % Show beam specific rows
      switch app.BeamTypeDropDown.Value
        case 'Gaussian'
          % Nothing to do
        case 'Hermite'
          app.ExtraGrid.RowHeight{app.MmodeSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.NmodeSpinner.Layout.Row} = rheight;
          app.MmodeSpinner.Enable = 'on';
          app.NmodeSpinner.Enable = 'on';
        case 'Laguerre'
          app.ExtraGrid.RowHeight{app.LmodeSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.PmodeSpinner.Layout.Row} = rheight;
          app.LmodeSpinner.Enable = 'on';
          app.PmodeSpinner.Enable = 'on';
        otherwise
          error('Internal error');
      end
    end
    
    function typeDropDownValueChanged(app)
      
      % Update visible shape widgets
      app.configureRowVisibility();
      
      % Continue processing changed value
      app.updateParametersCb();
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);
      
      % Create grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 9);
      app.ExtraGrid.RowHeight{end} = '1x';
      app.ExtraGrid.RowSpacing = 1;
      
      % TODO: Hide unneeded components from base
        
      % Beam type dropdown
      app.BeamTypeDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
        'label', 'Beam Type');
      app.BeamTypeDropDown.Items = {'Gaussian', 'Laguerre', 'Hermite'};
      app.BeamTypeDropDown.Layout.Row = 1;
      app.BeamTypeDropDown.Layout.Column = 1;
      app.BeamTypeDropDown.ValueChangedFcn = @(~,~) app.typeDropDownValueChanged();
        
      % Lmode spinner
      app.LmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'L mode');
      app.LmodeSpinner.Layout.Row = 2;
      app.LmodeSpinner.Layout.Column = 1;
      app.LmodeSpinner.Step = 1;
      app.LmodeSpinner.Limits = [-Inf,Inf];
      app.LmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Pmode spinner
      app.PmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'P mode');
      app.PmodeSpinner.Layout.Row = 3;
      app.PmodeSpinner.Layout.Column = 1;
      app.PmodeSpinner.Step = 1;
      app.PmodeSpinner.Limits = [0,Inf];
      app.PmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Mmode spinner
      app.MmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'M mode');
      app.MmodeSpinner.Layout.Row = 4;
      app.MmodeSpinner.Layout.Column = 1;
      app.MmodeSpinner.Step = 1;
      app.MmodeSpinner.Limits = [0,Inf];
      app.MmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Nmode spinner
      app.NmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'N mode');
      app.NmodeSpinner.Layout.Row = 5;
      app.NmodeSpinner.Layout.Column = 1;
      app.NmodeSpinner.Step = 1;
      app.NmodeSpinner.Limits = [0,Inf];
      app.NmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
        
      % NA spinner
      app.NaSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Numerical Aperture');
      app.NaSpinner.Layout.Row = 6;
      app.NaSpinner.Layout.Column = 1;
      app.NaSpinner.Step = 0.1;
      app.NaSpinner.LowerLimitInclusive = 'off';
      app.NaSpinner.Limits = [0.0, Inf];
      app.NaSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
        
      % Polarisation jones vector entry
      app.PolarisationEntry = ott.ui.support.JonesPolarisationEntry(...
          app.ExtraGrid);
      app.PolarisationEntry.Layout.Row = 7;
      app.PolarisationEntry.Layout.Column = 1;
      app.PolarisationEntry.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Power spinner  
      app.PowerSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Power');
      app.PowerSpinner.Layout.Row = 8;
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