classdef PowerSpectrum < matlab.apps.AppBase ...
    & ott.ui.support.AppProperties

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
  end
  
end