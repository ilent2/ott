classdef PowerSpectrum < ott.ui.support.AppTopLevel ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Generate plot of power spectrum from position time series data
%
% This app cab be launched from the app launcher via Tools -> PowerSpectrum
% or instantiated on the command line with:
%
%   ott.ui.tools.PowerSpectrum()

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Add support for input arguments to class

  properties (Constant)
    cnameText = 'PowerSpectrum';

    nameText = 'Generate Power Spectrum';

    aboutText = ['Generate plots of power spectrum from position time' ...
      'series data.'];
    
    helpText = {ott.ui.tools.PowerSpectrum.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 420];
  end
  
  properties
    MainGrid          matlab.ui.container.GridLayout
    XInputDropdown    ott.ui.support.VariableDropdown
    YInputDropdown    ott.ui.support.VariableDropdown
    XLabelEntry       ott.ui.support.LabeledTextEntry
    YLebelEntry       ott.ui.support.LabeledTextEntry
    UpdateButton      ott.ui.support.UpdateCheckButton
    UIAxes            matlab.ui.control.UIAxes
  end
  
  methods (Access=protected)
    function createMainComponents(app)
      
      % Main layout grid
      app.MainGrid = uigridlayout(app.UIFigure);
      app.MainGrid.RowHeight = repmat({32}, 1, 6);
      app.MainGrid.RowHeight{end-1} = '1x';
      app.MainGrid.ColumnWidth = {200, '1x'};
      app.MainGrid.ColumnSpacing = 30;
      app.MainGrid.RowSpacing = 0;
      
      % TODO: Rest of inputs
      
    end
  end
  
  methods (Access=public)
    function app=PowerSpectrum(t, x)
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTopLevel();
      
      if nargin ~= 0
        % TODO
      end
      
      if nargout == 0
        clear app;
      end
    end
  end
end