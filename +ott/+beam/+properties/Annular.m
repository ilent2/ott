classdef Annular < ott.beam.properties.MaterialBeam
% Properties describing a far-field annular beam.
% Inherits from :class:`MaterialBeam`.
%
% Properties
%   - angles      -- Angles defining range of annular (radians)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    angles        % Angles defining range of annular (radians).
  end

  methods % Getters/setters
    function beam = set.angles(beam, val)
      assert(isnumeric(val) && numel(val) == 2 ...
          && min(val) >= 0.0 && max(val) <= pi, ...
          'angles must be 2 element numeric in range [0, pi]');
      beam.angles = val;
    end
  end
end
