classdef Dda < ott.ui.tmatrix.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Dda';

    nameText = 'DDA T-matrix';

    aboutText = ['Generate T-matrix using the discrete dipole ' ...
      'approximation.'];
    
    helpText = {ott.ui.tmatrix.Dda.aboutText, ...
      ''};
  end
  
end