classdef (Abstract) AppTwoColumn < ott.ui.support.AppTopLevel
% Creates a two column application

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  % Properties that correspond to app components
  properties (Access = public)
    GridLayout                      matlab.ui.container.GridLayout
    LeftPanel                       matlab.ui.container.Panel
    RightPanel                      matlab.ui.container.Panel
  end
  
  % Properties that correspond to apps with auto-reflow
  properties (Access = private)
    onePanelWidth = 576;
    onePanelHeight = 300;
  end

  % Callbacks that handle component events
  methods (Access = private)

    % Changes arrangement of the app based on UIFigure width
    function updateAppLayout(app, ~)
      currentFigureWidth = app.UIFigure.Position(3);
      currentFigureHeight = app.UIFigure.Position(4);
      if(currentFigureWidth <= app.onePanelWidth ...
          || currentFigureHeight <= app.onePanelHeight)
        % Change to a 2x1 grid
        app.GridLayout.RowHeight = {418, 418};
        app.GridLayout.ColumnWidth = {'1x'};
        app.RightPanel.Layout.Row = 2;
        app.RightPanel.Layout.Column = 1;
      else
        % Change to a 1x2 grid
        app.GridLayout.RowHeight = {'1x'};
        app.GridLayout.ColumnWidth = {251, '1x'};
        app.RightPanel.Layout.Row = 1;
        app.RightPanel.Layout.Column = 2;
      end
    end
  end
  
  methods (Abstract, Access=protected)
    createLeftComponents(app)
    createRightComponents(app)
  end
  
  methods (Access=protected)
    % Create UIFigure and components
    function createMainComponents(app)
      
      % Configure figure
      app.UIFigure.AutoResizeChildren = 'off';
      app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);
      
      % Create GridLayout
      app.GridLayout = uigridlayout(app.UIFigure);
      app.GridLayout.ColumnWidth = {251, '1x'};
      app.GridLayout.RowHeight = {'1x'};
      app.GridLayout.ColumnSpacing = 0;
      app.GridLayout.RowSpacing = 0;
      app.GridLayout.Padding = [0 0 0 0];
      app.GridLayout.Scrollable = 'on';
      
      % Create LeftPanel
      app.LeftPanel = uipanel(app.GridLayout);
      app.LeftPanel.Layout.Row = 1;
      app.LeftPanel.Layout.Column = 1;
      
      app.createLeftComponents();
      
      % Create RightPanel
      app.RightPanel = uipanel(app.GridLayout);
      app.RightPanel.Layout.Row = 1;
      app.RightPanel.Layout.Column = 2;
      
      app.createRightComponents();
      
    end
  end
end