classdef Simple < ott.ui.shape.AppBase
% Generate a simple geometric shape.
% For more complex shapes, consider using an external CAD program
% and :class:`CadFileLoader`.
%
% Supported shapes:
%   - Sphere
%   - Ellipsoid
%   - Cylinder
%   - Cube
%   - Rectangular Prism
%   - Pill
%   - Bicone
%   - Cone Tipped Cylinder
%   - Biconcave Disc
%
% Some of these shapes may move to their own GUI in a future release.
%
% This GUI can be launched from the launcher under
% Shape -> Simple or running the following command:
%
%   ott.ui.shape.Simple()
%
% The shape is stored internally and/or written to the matlab workspace
% if a variable name is given for the shape.  To access the internal
% instance use:
%
%   app = ott.ui.shape.Simple()
%   shape = app.shape

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Constant)
    cnameText = 'Simple';

    nameText = 'Simple Shape Generator';

    aboutText = ['Generates a simple geometric shape.  For more complex' ...
      ' shapes consider using a external CAD program and loading the' ...
      ' shape using CadFileLoader.'];
    
    helpText = {ott.ui.shape.Simple.aboutText, ...
      ''};
  end
  
end