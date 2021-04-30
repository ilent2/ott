classdef Pointmatch < ott.ui.tmatrix.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Pointmatch';

    nameText = 'Pointmatch T-matrix';

    aboutText = ['Generate T-matrix for star shaped particle using' ...
      'the point matching method.'];
    
    helpText = {ott.ui.tmatrix.Pointmatch.aboutText, ...
      ''};
  end
  
end