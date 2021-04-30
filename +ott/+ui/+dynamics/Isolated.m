classdef Isolated < matlab.apps.AppBase ...
    & ott.ui.support.AppProperties

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Isolated';

    nameText = 'Simulate Isolated Dynamics';

    aboutText = ['Simulate dynamics of an isolated particle.'];
    
    helpText = {ott.ui.dynamics.Isolated.aboutText, ...
      ''};
  end
  
end