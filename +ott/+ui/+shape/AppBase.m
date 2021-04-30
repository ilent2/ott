classdef (Abstract) AppBase < matlab.apps.AppBase ...
    & ott.ui.support.AppProperties
% Base class for shape creation application windows.
%
% This class is not intended to be instantiated directly.
% The shape is stored internally and/or written to the matlab workspace
% if a variable name is given for the shape.  To access the internal
% instance use:
%
%   app = ott.ui.shape.<name-of-your-app>()
%   shape = app.shape

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    shape       % Internal representation of the shape
  end

end