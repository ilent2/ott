classdef Visualise < ott.ui.beam.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Visualise';

    nameText = 'Generate Beam Visualisation';

    aboutText = ['Generate visualisation of a beam.'];
    
    helpText = {ott.ui.beam.Visualise.aboutText, ...
      ''};
  end
  
end