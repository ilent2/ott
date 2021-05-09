classdef Mie < ott.ui.tmatrix.NewTmatrixBase
% Generate a T-matrix cotaining the Mie coefficients.

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
    
    windowName = ott.ui.tmatrix.Mie.nameText;
  end
  
  properties (Access=public)
    
  end
  
  methods (Access=protected)
    function data = generateData(app)
      data = ott.tmatrix.Mie(radius);  % TODO: parameters
    end
    
    function code = generateCode(app)
      code = {}; % TODO
    end
    
    function createMainComponents(app)
      % Create components for gui
      
      % Create base components
      createMainComponents@ott.ui.tmatrix.NewTmatrixBase(app);
      
      % TODO: Configure extra grid
      
      % Properties selector
%       app.PropertiesDropDown = ott.ui.support.LabeledDropDown(app.ExtraGrid, ...
%         'label', 'Properties');
%       % TODO: Layout
%       
%       % Radius entry
%       app.RadiusSpinner = ott.ui.support.LabeledSpinner(app.ExtraGrid, ...
%         'label', 'Radius');
%       % TODO: Layout
      
    end
  end
  
  methods (Access=public)
    function app=Mie()
      % Start the ForcePosition GUI
      
      app = app@ott.ui.tmatrix.NewTmatrixBase();
      
      if nargout == 0
        clear app;
      end
    end
  end
  
end