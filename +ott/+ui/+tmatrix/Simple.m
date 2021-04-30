classdef Simple < ott.ui.tmatrix.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Generate Simple T-matrix';

    aboutText = ['Generate T-matrix using a guess at an appropriate ' ...
      'method.  For spheres, defaults to Mie.  Other shapes use ' ...
      'point matching or DDA depending on the shape properties.'];
    
    helpText = {ott.ui.tmatrix.Simple.aboutText, ...
      ''};
  end
  
end