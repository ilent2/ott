classdef Dda < ott.ui.tmatrix.NewTmatrixBase
% Generate a T-matrix using the discrete dipole approximation.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Dda';

    nameText = 'DDA T-matrix';

    aboutText = ['Generate a T-matrix using the discrete dipole ' ...
      'approximation.'];
    
    helpText = {ott.ui.tmatrix.Dda.aboutText, ...
      '', ...
      'Spacing : Distance between voxels.', ...
      '', ...
      'Low Memory : When true, uses the low memory DDA method.'};
    
    windowName = ott.ui.tmatrix.Dda.nameText;
  end
  
  properties (Access=public)
    SpacingSpinner      ott.ui.support.LabeledSpinner
    LowMemoryCheckhbox  matlab.ui.control.CheckBox
  end
  
  methods (Access=protected)
    
    function setDefaultValues(app)
      % Set default values of our widgets
      
      app.SpacingSpinner.Value = 1/20;
      app.LowMemoryCheckhbox.Value = false;
      
      % Let base do the rest
      setDefaultValues@ott.ui.tmatrix.NewTmatrixBase(app);
    end
    
    function data = generateData(app)
      % Generate T-matrix
      
      data = ott.tmatrix.Dda.FromShape(app.ShapeName.Variable, ...
          'spacing', app.SpacingSpinner.Value, ...
          'index_relative', app.RelativeIndexSpinner.Value, ...
          'low_memory', app.LowMemoryCheckbox.Value);
    end
    
    function code = generateCode(app)
      code = {};
      
      % Code for T-matrix generation
      code{end+1} = ['tmatrix = ott.tmatrix.Dda.FromShape(' app.ShapeName.Value ', ...'];
      code{end+1} = ['    ''spacing'', ' app.SpacingSpinner.Value ', ...'];
      code{end+1} = ['    ''index_relative'', ' app.RelativeIndexSpinner.Value ', ...'];
      code{end+1} = ['    ''low_memory'', ' app.LowMemoryCheckbox.Value ');'];
    end
    
    function createMainComponents(app)
      
      % Call base for most things
      createMainComponents@ott.ui.tmatrix.NewTmatrixBase(app);
      
      % Create grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 3);
      app.ExtraGrid.RowHeight{end} = '1x';
      app.ExtraGrid.RowSpacing = 1;
      
      % Spacing spinner
      app.SpacingSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'label', 'Spacing');
      app.SpacingSpinner.Layout.Row = 1;
      app.SpacingSpinner.Layout.Column = 1;
      app.SpacingSpinner.Step = 0.01;
      app.SpacingSpinner.Limits = [0,Inf];
      app.SpacingSpinner.LowerLimitInclusive = false;
      app.SpacingSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
      % Low memory checkbox
      app.LowMemoryCheckhbox = uicheckbox(app.ExtraGrid);
      app.LowMemoryCheckhbox.Text = 'Low Memory';
      app.LowMemoryCheckhbox.Layout.Column = 1;
      app.LowMemoryCheckhbox.Layout.Row = 2;
      app.LowMemoryCheckhbox.ValueChangedFcn = @(~,~) app.updateParametersCb();
      
    end
  end
  
  methods (Access=public)
    function app=Dda()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.NewTmatrixBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end