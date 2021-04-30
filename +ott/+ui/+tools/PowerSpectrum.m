classdef PowerSpectrum < ott.ui.support.AppTopLevel
% Generate plot of power spectrum from position time series data
%
% This app cab be launched from the app launcher via Tools -> PowerSpectrum
% or instantiated on the command line with:
%
%   ott.ui.tools.PowerSpectrum
%
% The constructor also supports an optional input for the data to
% calculate the time series from:
%
%   ott.ui.tools.PowerSpectrum(t, x)

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

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
  
  methods (Access=protected)
    function startupFcn(app)
    end
    
    function createComponents(app)
      createComponents@ott.ui.support.AppTopLevel(app);
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