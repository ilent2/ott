classdef Simple < ott.ui.support.AppTopLevel

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Generate Simple Particle';

    aboutText = ['Generate a simple particle representation.'];
    
    helpText = {ott.ui.particle.Simple.aboutText, ...
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
    function app=Simple()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTopLevel();
      if nargout == 0
        clear app;
      end
    end
  end
  
end