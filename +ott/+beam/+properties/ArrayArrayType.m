classdef ArrayArrayType < ott.beam.properties.ArrayType
% Declares a constant array_type property for array array types.
%
% This class is intended for use in ott.Beam instances which only
% support a single array type.  For beams which support multiple
% array types, use :class:`AnyArrayType` instead.
%
% Properties
%   - array_type      -- Constant, 'array'

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    array_type = 'array';
  end

  methods % Getters/setters
    function beam = set.array_type(beam, val)
      error('Array type can not be changed for this beam');
    end
  end
end
