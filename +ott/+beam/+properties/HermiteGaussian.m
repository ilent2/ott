classdef Scattered < ott.beam.properties.Beam
% Properties of scattered beams
%
% Properties
%   - waist      -- Beam waist at focus
%   - mmode      -- Hermite mode order
%   - nmode      -- Hermite mode order
%   - power         -- Beam power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    mmode       % Hermite mode order
    nmode       % Hermite mode order
  end

  methods % Getters/setters
    % mmode       % Hermite mode order
    % nmode       % Hermite mode order

    function beam = set.nmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'nmode must be numeric scalar');
      assert(round(val) == val, 'nmode must be integer');
      beam.nmode = val;
    end

    function beam = set.mmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'mmode must be numeric scalar');
      assert(round(val) == val, 'mmode must be integer');
      beam.mmode = val;
    end
  end
end
