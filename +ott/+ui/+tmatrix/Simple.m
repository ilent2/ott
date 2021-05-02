classdef Simple < ott.ui.tmatrix.NewTmatrixBase
% Generate a simple T-matrix by taking a guess at the appropriate
% method based on the shape geometry.

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Generate Simple T-matrix';

    aboutText = ['Generate T-matrix using a guess at an appropriate ' ...
      'method.  For spheres, defaults to Mie.  Other shapes use ' ...
      'point matching or DDA depending on the shape properties.'];
    
    helpText = {ott.ui.tmatrix.Simple.aboutText, ...
      ''};
    
    windowName = ott.ui.tmatrix.Simple.nameText;
  end
  
  methods (Access=protected)
    function startupFcn(app)
    end
  end
  
  methods (Access=public)
    function app=Simple()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.NewTmatrixBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
end