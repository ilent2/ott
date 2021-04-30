classdef Simple < matlab.apps.AppBase ...
    & ott.ui.support.AppProperties

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Generate Simple Particle';

    aboutText = ['Generate a simple particle representation.'];
    
    helpText = {ott.ui.particle.Simple.aboutText, ...
      ''};
  end
  
end