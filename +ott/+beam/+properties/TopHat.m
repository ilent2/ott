classdef TopHat
% Properties of a collimated Top-Hat beam.
%
% Properties
%   - radius        -- Radius of the top-hat

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    radius
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.TopHat')
        args = ott.utils.addDefaultParameter('radius', other.radius, args);
      end
    end
  end

  methods
    function beam = TopHat(varargin)
      % Construct TopHat properties
      %
      % Usage
      %   beam = beam@TopHat(radius)

      p = inputParser;
      p.addOptional('radius', [], @isnumeric);
      p.parse(varargin{:});

      beam.radius = p.Results.radius;
    end
  end

  methods % Getters/setters
    function beam = set.radius(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'radius must be numeric scalar');
    end
  end
end
