classdef BscPointmatch < ott.Bsc
%BscPointmatch base class for BSC generated using point matching
% Provides support for both farfield and focal plane point matching.
%
% BscPointmatch properties:
%
% BscPointmatch methods:
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function beam = BscPointmatch(varargin)
      beam = beam@ott.Bsc(varargin{:});
    end
  end
end
