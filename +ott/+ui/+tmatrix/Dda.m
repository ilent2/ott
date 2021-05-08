classdef Dda < ott.ui.tmatrix.NewTmatrixBase
% Generate a T-matrix using the discrete dipole approximation.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Dda';

    nameText = 'DDA T-matrix';

    aboutText = ['Generate a T-matrix using the discrete dipole ' ...
      'approximation.'];
    
    helpText = {ott.ui.tmatrix.Dda.aboutText, ...
      ''};
    
    windowName = ott.ui.tmatrix.Dda.nameText;
  end
  
  methods (Access=protected)
    function data = GenerateData(app)
      data = ott.tmatrix.Dda();  % TODO: parameters
    end
  end
  
  methods (Access=public)
    function app=Dda()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.NewTmatrixBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end