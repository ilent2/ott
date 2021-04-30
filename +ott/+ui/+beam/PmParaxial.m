classdef PmParaxial < ott.ui.beam.AppBase

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'PmParaxial';

    nameText = 'Paraxial Pointmatched Beam';

    aboutText = ['Generates a beam using paraxial point matching from' ...
      ' a given E field profile.'];
    
    helpText = {ott.ui.beam.PmParaxial.aboutText, ...
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
    function app=PmParaxial()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.beam.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end