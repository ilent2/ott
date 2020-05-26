classdef AnyArrayType < ott.beam.properties.ArrayType
% Declares a variable array_type property for beams supporting multiple types.
%
% This class by default supports coherent, incoherent or array values
% for the array_type property.  It is intended for use by beams supporting
% multiple array types.  The options for array_type can be further restricted
% by overloading the 'validateArrayType' property.
%
% Methods
%   - validateArrayType   -- Validates the array type
%
% Properties
%   - array_type      -- Array type, 'coherent', 'incoherent' or 'array'

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    array_type
  end

  methods
    function beam = AnyArrayType(varargin)
      % Store array type argument from inputs
      %
      % Usage
      %   beam = beam@ott.beam.utils.AnyArray(array_type)
      %   Can also use named arguments.

      p = inputParser;
      p.addOptional('array_type', [], ...
          @(x) any(strcmpi('coherent', 'array', 'incoherent')));
      p.parse(varargin{:});
      beam.array_type = p.Results.array_type;
    end
  end

  methods (Hidden)
    function validateArrayType(newtype)
      % Method called to validate type changes
      %
      % Default method has not validation.
      % Overload this method to limit type changes.
    end
  end

  methods % Getters/setters
    function beam = set.array_type(beam, val)

      % Check valid value
      assert(any(strcmpi(val, {'coherent', 'array', 'incoherent'})), ...
        'array_type must be one of ''coherent'', ''array'', ''incoherent''');

      % Apply additional validations
      beam.validateArrayType();

      beam.array_type = val;
    end
  end
end

