classdef AppProperties < handle
% Defines the properties expected by apps used in the Launcher.
%
% Abstract properties:
%   - cnameText
%   - nameText
%   - aboutText
%   - helpText

% Copyright 2021 IST Austria, Written by Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract, Constant)
    cnameText
    nameText
    aboutText
    helpText
  end
end
