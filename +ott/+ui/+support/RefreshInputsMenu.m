classdef RefreshInputsMenu < handle
% Adds a Refresh Inputs menu item to the File menu\
%
% Methods
%   - registerInput

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Access=protected)
    inputs (1, :) ott.ui.support.VariableDropDown
  end
  
  methods (Access=protected)
    function registerRefreshInput(app, input)
      % Store reference to input
      app.inputs(end+1) = input;
    end
  end
  
  methods (Access={?ott.ui.support.AppTopLevel})
    function refreshInputsCb(app, ~)
      
      for inp = app.inputs
        inp.refresh();
      end
    end
  end
end