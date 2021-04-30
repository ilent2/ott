classdef Scattered < ott.ui.beam.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Scattered';

    nameText = 'Scatter Beam';

    aboutText = ['Generate a scattered beam from a particle and incident' ...
      'beam instance.'];
    
    helpText = {ott.ui.beam.Scattered.aboutText, ...
      ''};
  end
  
end