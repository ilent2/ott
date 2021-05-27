classdef Pointmatch < ott.ui.tmatrix.NewTmatrixBase
% Generate a T-matrix using the point matching method.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Pointmatch';

    nameText = 'Pointmatch T-matrix';

    aboutText = ['Generate T-matrix for star shaped particle using' ...
      'the point matching method.'];
    
    helpText = {ott.ui.tmatrix.Pointmatch.aboutText, ...
      ''};
    
    windowName = ott.ui.tmatrix.Pointmatch.nameText;
  end
  
  properties
    NmaxSpinner       ott.ui.support.LabeledSpinner
  end
  
  methods (Access=protected)
    
    function setDefaultValues(app)
      % Set default values of our widgets
      
      app.NmaxSpinner.Value = 20;
      
      % Let base do the rest
      setDefaultValues@ott.ui.tmatrix.NewTmatrixBase(app);
    end
    
    function data = generateData(app)
      data = ott.tmatrix.Pointmatch.FromShape(...
        app.ShapeName.Variable ./ app.WavelengthSpinner.Value, ...
        app.RelativeIndexSpinner.Value, ...
        'Nmax', app.NmaxSpinner.Value);
    end
    
    function code = generateCode(app)
      % Generate code for app
      
      code = {};
      
      code{end+1} = ['Nmax = ' num2str(app.NmaxSpinner.Value) ';'];
      code{end+1} = ['wavelength = ' num2str(app.WavelengthSpinner.Value) ';'];
      code{end+1} = ['index_relative = ' num2str(app.RelativeIndexSpinner.Value) ';'];
      code{end+1} = 'tmatrix = ott.tmatrix.Pointmatch.FromShape(...';
      code{end+1} = '  shape./wavelength, index_relative, ''Nmax'', Nmax);';
    end
    
    function shapeChangedCb(app)
      % Update the Nmax spinner to a new default when shape changes.
      try
        shape = app.ShapeName.Variable;
        wavelength = app.WavelengthSpinner.Value;
        Nmax = ott.utils.ka2nmax(2*pi*shape.maxRadius./wavelength);
        app.NmaxSpinner.Value = Nmax;
      catch ME
        warning('ott:ui:pointmatch:shape_update_nmax', ...
          ['Couldn''t set new Nmax from shape variable:' newline, ...
          ME.message]);
      end
    end
    
    function shapeWavelengthChangedCb(app)
      app.shapeChangedCb();
      app.updateParametersCb();
    end
    
    function createMainComponents(app)
      
      % Call base for most things
      createMainComponents@ott.ui.tmatrix.NewTmatrixBase(app);
      
      % Create grid
      app.ExtraGrid.RowHeight = repmat({32}, 1, 2);
      app.ExtraGrid.RowHeight{end} = '1x';
      app.ExtraGrid.RowSpacing = 1;
      
      % Add callback for shape name change and wavelength changed
      app.ShapeName.ValueChangedFcn = @(~,~) app.shapeChangedCb();
      app.WavelengthSpinner.ValueChangedFcn = @(~,~) app.shapeWavelengthChangedCb();
      
      % Nmax spinner
      app.NmaxSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
          'wwidth', app.wwidth, 'label', 'Nmax');
      app.NmaxSpinner.Layout.Row = 1;
      app.NmaxSpinner.Layout.Column = 1;
      app.NmaxSpinner.Step = 1;
      app.NmaxSpinner.Limits = [0,Inf];
      app.NmaxSpinner.LowerLimitInclusive = false;
      app.NmaxSpinner.ValueChangedFcn = @(~,~) app.updateParametersCb();
    end
  end
  
  methods (Access=public)
    function app=Pointmatch()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.NewTmatrixBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end