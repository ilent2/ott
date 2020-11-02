classdef Mapping
% Declares a mapping property for paraxial to far-field mapped beams
%
% Properties
%   - mapping       -- Paraxial to far-field beam mapping

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    mapping  % Paraxial to far-field beam mapping
  end

  properties (Hidden, SetAccess=protected)
    mappingInternal
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.mapping(beam, val)
      assert(any(strcmpi(val, {'sin', 'tan', 'theta'})), ...
          'mapping must be one of ''sin'' ''tan'' ''theta''');
      beam.mappingInternal = val;
      beam.data = [];
    end
    function val = get.mapping(beam)
      val = beam.mappingInternal;
    end
  end
end

