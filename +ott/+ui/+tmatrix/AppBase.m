classdef AppBase < ott.ui.support.AppTopLevel

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    tmatrix
  end
  
  properties (Access=public)
    MainGrid            matlab.ui.container.GridLayout
    VariableName        ott.ui.support.OutputVariableEntry
    ParticleName        ott.ui.support.VariableDropdown
    UpdateButton        ott.ui.support.UpdateWithProgress
  end
  
  methods (Access=protected)
    function createMainComponents(app)
      
      % Create grid
      app.MainGrid = uigridlayout(app.UIFigure, [4, 1]);
      app.MainGrid.RowHeight = {32, 32, '1x', 32};
      app.MainGrid.ColumnWidth = {220};
      app.MainGrid.ColumnSpacing = 1;
      app.MainGrid.RowSpacing = 1;
      
      app.VariableName = ott.ui.support.OutputVariableEntry(app.MainGrid);
      app.VariableName.Layout.Row = 1;
      app.VariableName.Layout.Column = 1;
      
      app.ParticleName = ott.ui.support.VariableDropdown(app.MainGrid);
      app.ParticleName.Layout.Row = 2;
      app.ParticleName.Layout.Column = 1;
      
      app.UpdateButton = ott.ui.support.UpdateWithProgress(app.MainGrid);
      app.UpdateButton.Layout.Row = 4;
      app.UpdateButton.Layout.Column = 1;
      
    end
  end
  
end