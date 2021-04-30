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
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 420];
  end
  
  methods (Access=protected)
    function startupFcn(app)
    end
    
    function createComponents(app)
      createComponents@ott.ui.support.AppTopLevel(app);
    end
  end
  
  methods (Access=public)
    function app=Dda()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end