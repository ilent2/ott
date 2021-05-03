classdef (Abstract) NewTmatrixBase < ott.ui.support.AppTopLevel
% Base class for T-matrix creation applications.
%
% Properties:
%   - tmatrix

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    windowSize = [640, 420];
  end

  properties (SetAccess=protected)
    tmatrix
  end
  
  properties (Access=public)
    MainGrid            matlab.ui.container.GridLayout
    VariableName        ott.ui.support.OutputVariableEntry
    ParticleName        ott.ui.support.VariableDropdown
    UpdateButton        ott.ui.support.UpdateWithProgress
  end
  
  methods (Access=protected, Abstract)
    generateTmatrix(app)
    generateCode(app)
  end
  
  methods (Access=protected)
    function setDefaultValues(app)
      app.VariableName.Value = '';
      app.ParticleName.Value = '';
      app.UpdateButton.Value = 0;
      app.UpdateButton.ClearErrors();
    end
    
    function updateCb(app, ~)
      % Called when a value is changed or when update is clicked
      
      % Generate new beam
      app.tmatrix = app.generateTmatrix();
      
      % Write to workspace (T-matrix doesn't support preview)
      app.VariableName.WriteVariable(app.beam);
      
    end
    
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