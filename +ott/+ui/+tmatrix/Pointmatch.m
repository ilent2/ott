classdef Pointmatch < ott.ui.tmatrix.NewTmatrixBase
% Generate a T-matrix using the point matching method.

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
    
    windowName = ott.ui.tmatrix.Pointmatch.nameText;
  end
  
  methods (Access=protected)
    function startupFcn(app)
    end
  end
  
  methods (Access=public)
    function app=Pointmatch()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.NewTmatrixBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end