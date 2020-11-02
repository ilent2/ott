classdef Annular
% Properties of annular beams
%
% Properties
%   - theta       -- Annular angle [radians]

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    theta
  end

  properties (Hidden, SetAccess=protected)
    thetaInternal
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.theta(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'theta must be numeric scalar');
      beam.thetaInternal = val;
      beam.data = [];
    end
    function val = get.theta(beam)
      val = beam.thetaInternal;
    end
  end
end
