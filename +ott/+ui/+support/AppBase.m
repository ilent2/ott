classdef (Abstract) AppBase < matlab.apps.AppBase
% Base class for all optical tweezers toolbox GUI apps.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Access=public)
    UIFigure                        matlab.ui.Figure
  end
  
  properties (Abstract, Constant)
    windowSize
    windowName
  end
  
  properties (Constant)
    windowOffset = [100, 100];
  end
  
  methods (Abstract, Access=protected)
    createComponents(app)
  end
  
  methods (Access=protected)
    function startupFcn(app)
      % Overload this method for startup functionality
    end
  end
  
  methods (Access=public)
    function app = AppBase()
      % Start the CadFileLoader interface
      
      % Create window, hide until after components created
      app.UIFigure = uifigure('Visible', 'off');
      app.UIFigure.Position = [app.windowOffset, app.windowSize];
      app.UIFigure.Name = app.windowName;

      % Create UI
      app.createComponents();
      
      % Make window visible
      app.UIFigure.Visible = 'on';
      
      % Execute the startup function
      runStartupFcn(app, @startupFcn)
      
      % Register the app with App Designer
      registerApp(app, app.UIFigure);

      if nargout == 0
        clear app;
      end
    end
    
    % Code that executes before app deletion
    function delete(app)

      % Delete UIFigure when app is deleted
      delete(app.UIFigure)
    end
  end
end