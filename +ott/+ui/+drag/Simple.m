classdef Simple < matlab.apps.AppBase ...
    & ott.ui.support.AppProperties

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Simple Drag';

    aboutText = ['Calculate drag tensor for a simple shape.'];
    
    helpText = {ott.ui.drag.Simple.aboutText, ...
      ''};
  end
  
end