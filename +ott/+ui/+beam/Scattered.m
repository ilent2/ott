classdef Scattered < ott.ui.beam.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Scattered';

    nameText = 'Scatter Beam';

    aboutText = ['Generate a scattered beam from a particle and incident' ...
      'beam instance.'];
    
    helpText = {ott.ui.beam.Scattered.aboutText, ...
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
    function app=Scattered()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end