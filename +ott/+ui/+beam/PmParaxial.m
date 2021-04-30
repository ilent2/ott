classdef PmParaxial < ott.ui.beam.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'PmParaxial';

    nameText = 'Paraxial Pointmatched Beam';

    aboutText = ['Generates a beam using paraxial point matching from' ...
      ' a given E field profile.'];
    
    helpText = {ott.ui.beam.PmParaxial.aboutText, ...
      ''};
  end
  
end