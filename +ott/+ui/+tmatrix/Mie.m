classdef Mie < ott.ui.tmatrix.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Mie';

    nameText = 'Mie T-matrix';

    aboutText = ['Generate T-matrix for a spherical particle using' ...
      'the Mie coefficients.'];
    
    helpText = {ott.ui.tmatrix.Mie.aboutText, ...
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
    function app=Mie()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end