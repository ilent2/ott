classdef BesselLike < ott.ui.beam.NewBeamBase ...
    & ott.ui.support.GenerateCodeMenu
% Generate a simple Bessel-like beam (Mathieu, Webber or Bessel).
%
% This GUI can be launched from the launcher Beam -> Simple or by running:
%
%   ott.ui.beam.BesselLike

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'BesselLike';

    nameText = 'Generate Bessel-like Beam';

    aboutText = ['Generate a bessel-like beam, such as a Mathieu, Webber', ...
      ' or Bessel beam.'];
    
    helpText = {ott.ui.beam.BesselLike.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.BesselLike.nameText;
  end
  
  properties (Access=public)
    TypeDropDown        ott.ui.support.LabeledDropDown
    ThetaSpinner        ott.ui.support.LabeledSpinner
    AlphaSpinner        ott.ui.support.LabeledSpinner
    EllipticitySpinner  ott.ui.support.LabeledSpinner
    LmodeSpinner        ott.ui.support.LabeledSpinner
    MorderSpinner       ott.ui.support.LabeledSpinner
    ParityDropDown      ott.ui.support.LabeledDropDown
    PolBasisDropDown    ott.ui.support.LabeledDropDown
    PolSpinners         ott.ui.support.JonesPolarisationEntry
  end
  
  methods (Access=protected)
    
    function code = generateCode(app)
      code = {}; % TODO
      
      switch app.TypeDropDown.Value
        case 'Bessel'
        case 'Webber'
        case 'Mathieu'
        otherwise
          error('Internal error');
      end
    end
    
    function data = generateData(app)
      data = []; % TODO
      
      switch app.TypeDropDown.Value
        case 'Bessel'
        case 'Webber'
        case 'Mathieu'
        otherwise
          error('Internal error');
      end
    end
    
    function setDefaultValues(app)
      
      app.TypeDropDown.Value = 'Bessel';
      app.ThetaSpinner.Value = 45;
      app.AlphaSpinner.Value = 1;
      app.EllipticitySpinner.Value = 1;
      app.MorderSpinner.Value = 0;
      app.LmodeSpinner.Value = 0;
      app.ParityDropDown.Value = 'Even';
      app.PolBasisDropDown.Value = 'Linear';
      app.PolSpinners.Value = [1, 1i];
      
      % Show/hide widgets
      app.configureRowVisibility();
      
      % Do base work
      setDefaultValues@ott.ui.beam.NewBeamBase(app);
    end
    
    function configureRowVisibility(app)
      % Update the visible widgets associated with the type selection
      
      rheight = 32;
      
      % Hide all widgets
      app.ExtraGrid.RowHeight(3:9) = repmat({0}, 1, 7);
      
      % Disable all widgets
      widgets = findobj(app.ExtraGrid.Children, '-property', 'Enable');
      set(widgets, 'Enable', 'off');
      app.TypeDropDown.Enable = 'on';
      app.ThetaSpinner.Enable = 'on';
      
      % Enable widgets for selected menu
      switch app.TypeDropDown.Value
        case 'Bessel'
          app.ExtraGrid.RowHeight{app.PolBasisDropDown.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.PolSpinners.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.LmodeSpinner.Layout.Row} = rheight;
          app.PolBasisDropDown.Enable = 'on';
          app.PolSpinners.Enable = 'on';
          app.LmodeSpinner.Enable = 'on';
        case 'Webber'
          app.ExtraGrid.RowHeight{app.AlphaSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.ParityDropDown.Layout.Row} = rheight;
          app.AlphaSpinner.Enable = 'on';
          app.ParityDropDown.Enable = 'on';
        case 'Mathieu'
          app.ExtraGrid.RowHeight{app.MorderSpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.EllipticitySpinner.Layout.Row} = rheight;
          app.ExtraGrid.RowHeight{app.ParityDropDown.Layout.Row} = rheight;
          app.MorderSpinner.Enable = 'on';
          app.EllipticitySpinner.Enable = 'on';
          app.ParityDropDown.Enable = 'on';
        otherwise
          error('Internal error');
      end
    end
    
    function typeChangedCb(app)
      
      % Show/hide widgets
      app.configureRowVisibility();
      
      % Continue to regular update stuff
      app.updateParametersCb();
    end
    
    function createLeftComponents(app)
      
      % Call base for most things
      createLeftComponents@ott.ui.beam.NewBeamBase(app);
      
      % Configure extra grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 10);
      app.ExtraGrid.RowHeight{end} = '1x';
      app.ExtraGrid.RowSpacing = 1;
      
      % TODO: Set base widget visibility
      
      % Type selector
      app.TypeDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
        'label', 'Type');
      app.TypeDropDown.Items = {'Bessel', 'Webber', 'Mathieu'};
      app.TypeDropDown.Layout.Row = 1;
      app.TypeDropDown.Layout.Column = 1;
      app.TypeDropDown.ValueChangedFcn = @(~,~) app.typeChangedCb();
      
      % Theta selector (all)
      app.ThetaSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Theta');
      app.ThetaSpinner.Layout.Row = 2;
      app.ThetaSpinner.Layout.Column = 1;
      app.ThetaSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Alpha spinner (webber)
      app.AlphaSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Alpha');
      app.AlphaSpinner.Layout.Row = 3;
      app.AlphaSpinner.Layout.Column = 1;
      app.AlphaSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Ellipticity spinner (Mathieu)
      app.EllipticitySpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Ellipticity');
      app.EllipticitySpinner.Layout.Row = 4;
      app.EllipticitySpinner.Layout.Column = 1;
      app.EllipticitySpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % L-mode spinner (Bessel)
      app.LmodeSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'Azim. Mode');
      app.LmodeSpinner.Layout.Row = 5;
      app.LmodeSpinner.Layout.Column = 1;
      app.LmodeSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % M-order spinner (Mathieu)
      app.MorderSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
        'label', 'M-order');
      app.MorderSpinner.Layout.Row = 6;
      app.MorderSpinner.Layout.Column = 1;
      app.MorderSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Parity (Mathieu/Webber)
      app.ParityDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
        'label', 'Parity');
      app.ParityDropDown.Items = {'Even', 'Odd'};
      app.ParityDropDown.Layout.Row = 7;
      app.ParityDropDown.Layout.Column = 1;
      app.ParityDropDown.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Polarisation basis drop down (Bessel)
      app.PolBasisDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
        'label', 'Pol. Basis');
      app.PolBasisDropDown.Items = {'Linear', 'Polar'};
      app.PolBasisDropDown.Layout.Row = 8;
      app.PolBasisDropDown.Layout.Column = 1;
      app.PolBasisDropDown.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Polarisation (Bessel)
      app.PolSpinners = ott.ui.support.JonesPolarisationEntry(app.ExtraGrid, ...
        'label', 'Polarisation');
      app.PolSpinners.Layout.Row = 9;
      app.PolSpinners.Layout.Column = 1;
      app.PolSpinners.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
    end
  end
  
  methods (Access=public)
    function app = BesselLike()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.NewBeamBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end