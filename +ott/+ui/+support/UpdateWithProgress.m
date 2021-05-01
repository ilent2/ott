classdef UpdateWithProgress < ott.ui.support.GridWidget
% Create an update button with an adjacent progress bar
  
% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    Gauge         matlab.ui.control.LinearGauge
    Button      matlab.ui.control.Button
  end
  
  methods
    function obj = UpdateWithProgress(parent, varargin)
      
      obj = obj@ott.ui.support.GridWidget(parent);
      
      % COnfigure grid
      obj.Grid.RowHeight = {22};
      obj.Grid.ColumnWidth = {70, 120};
      obj.Grid.ColumnSpacing = 10;
      obj.Grid.RowSpacing = 1;
      
      % Create BeamDropDownLabel
      obj.Button = uibutton(obj.Grid, 'push');
      obj.Button.Text = 'Update';
      obj.Button.Layout.Column = 1;
      obj.Button.Layout.Row = 1;

      % Create BeamDropDown
      obj.Gauge = uigauge(obj.Grid, 'linear');
      obj.Gauge.Enable = 'off';
      obj.Gauge.FontSize = 8;
      obj.Gauge.Layout.Column = 2;
      obj.Gauge.Layout.Row = 1;
      
    end
  end
end