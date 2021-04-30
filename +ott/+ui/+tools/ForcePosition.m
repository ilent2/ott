classdef ForcePosition < ott.ui.support.AppTopLevel

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'ForcePosition';

    nameText = 'Generate Force-Position plot';

    aboutText = ['Generate plots of the force as a function of ' ...
      'displacement or rotation angle.'];
    
    helpText = {ott.ui.tools.ForcePosition.aboutText, ...
      ''};
    
    windowName = ott.ui.tools.ForcePosition.nameText;
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
    function app=ForcePosition()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTopLevel();
      if nargout == 0
        clear app;
      end
    end
  end
end
