classdef Empty < ott.beam.properties.MaterialBeam & ott.beam.utils.ZeroPower
% Properties of scattered beams
%
% Properties
%   power         -- Constant, 0 (empty beams have no power)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Empty(varargin)
      beam = beam@ott.beam.properties.Beam(varargin{:});
    end
  end
end
