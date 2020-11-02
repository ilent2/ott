classdef Lmode
% Declares a lmode property for beams with angular momentum.
%
% Properties
%   - lmode       -- (numeric) Orbital angular momentum of the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    lmode       % Orbital angular momentum of the beam
  end

  properties (Hidden, SetAccess=protected)
    lmodeInternal;
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.lmode(beam, val)
      assert(isnumeric(val) && isscalar(val) && round(val) == val, ...
          'lmode must be numeric integer');
      beam.lmodeInternal = val;
      beam.data = [];
    end
    function val = get.lmode(beam)
      val = beam.lmodeInternal;
    end
  end
end
