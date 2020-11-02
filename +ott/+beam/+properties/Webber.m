classdef Webber < ott.beam.properties.Annular & ott.beam.properties.Parity
% Properties of a Webber beam
%
% Properties
%   - theta       -- Annular angle [radians]
%   - alpha       -- Parameter describe Webber beam
%   - parity      -- Parity of beam (either 'even' or 'odd')

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    alpha        % Parameter describe Webber beam
  end

  properties (Hidden, SetAccess=protected)
    alphaInternal
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.alpha(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'alpha must be numeric scalar');
      beam.alphaInternal = val;
      beam.data = [];
    end
    function val = get.alpha(beam)
      val = beam.alphaInternal;
    end
  end
end

