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
    function app=Pointmatch()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end