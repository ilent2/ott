classdef Visualise < ott.ui.support.AppTwoColumn ...
    & ott.ui.support.GenerateCodeMenu ...
    & ott.ui.support.RefreshInputsMenu
% Graphical user interface to visualise a beam

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Visualise';

    nameText = 'Generate Beam Visualisation';

    aboutText = ['Generate visualisation of a beam.'];
    
    helpText = {ott.ui.beam.Visualise.aboutText, ...
      ''};
    
    windowName = ott.ui.beam.PmParaxial.nameText;
    windowSize = [640, 420];
  end
  
  methods (Access=protected)
    
    function setDefaultValues(app)
      % TODO
    end
    
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function createRightComponents(app)
    end
    
    function createLeftComponents(app)
    end
  end
  
  methods (Access=public)
    function app=Visualise()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.support.AppTwoColumn();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end