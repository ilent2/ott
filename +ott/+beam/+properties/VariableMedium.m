classdef VariableMedium
% Adds a variable medium property to a Beam.
%
% Doesn't provide a constructor, uses properties.Beam.
%
% Properties
%   - medium
%
% Static methods
%   - likeProperties    -- Adds the medium to the properties list

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    medium
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      if isa(other, 'ott.beam.properties.VariableMedium')
        args = ott.utils.addDefaultParameter('medium', other.medium, args);
      end
    end
  end

  methods % Getters/setters
    function beam = set.medium(beam, val)
      assert(isa(val, 'ott.beam.medium.Medium') && numel(val) == 1, ...
          'medium must be a single ott.beam.medium.Medium');
      beam.medium = val;
    end
  end
end
