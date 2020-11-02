classdef Array < ott.beam.Beam
% Base class for coherent/incoherent arrays of beams.
% Inherits from :class:`Beam`.
%
% Properties
%   - data        -- Array of coherent/incoherent beams
%
% Abstract methods
%   - validateArray   -- Validate array data

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    data      % Array of coherent/incoherent beams
  end

  methods (Abstract, Hidden)
    validateArray
  end

  methods
    function beam = Array(data)
      % Construct instance of ott.beam.Array
      beam.data = data;
    end
  end

  methods % Getters/setters
    function beam = set.data(beam, val)
      assert(isa(val, 'ott.beam.Beam'), ...
          'data must be a array of ott.beam.Beam');
      beam.validateArray(data);
      beam.data = val;
    end
  end
end
