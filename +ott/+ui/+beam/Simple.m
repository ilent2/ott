classdef Simple < ott.ui.beam.AppBase
% Generate a simple beam representation and visualise.
%
% Supported beams:
%   - Gaussian
%   - Laguerre-Gaussian
%   - Hermite-Gaussian
%   - Plane Wave
%   - Mathieu
%   - Webber
%   - Bessel
%
% Some of these beams might move to their own interface in a future
% release.
%
% This GUI can be launched from the launcher Beam -> Simple or by running:
%
%   ott.ui.beam.Simple

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Should we split this interface into Gaussian, PlaneWave
% and Annular?

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Simple Beam';

    aboutText = ['Generate a simple beam.'];
    
    helpText = {ott.ui.beam.Simple.aboutText, ...
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
      
      app = app@ott.ui.beam.AppBase();
      if nargout == 0
        clear app;
      end
    end
  end
  
end