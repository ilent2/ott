classdef Bessel < ott.beam.properties.MaterialBeam
% Properties of a Bessel beam.
%
% Properties
%   - angle       - Far-field angle of Bessel beam (radians)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    angle       % Far-field angle of Bessel beam (radians)
  end

  methods % Getters/setters
    function beam = set.angle(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'angle must be numeric scalar');
      beam.angle = val;
    end
  end
end
