classdef Annular
% Properties describing a far-field annular beam.
%
% Properties
%   - angles      -- Angles defining range of annular (radians)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    angles        % Angles defining range of annular (radians).
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.Annular')
        args = ott.utils.addDefaultParameter('angles', other.angles, args);
      end
    end
  end

  methods
    function beam = Annular(varargin)
      % Construct Annular properties
      %
      % Usage
      %   beam = beam@FocussedTopHat(angles)

      p = inputParser;
      p.addOptional('angles', [], @isnumeric);
      p.parse(varargin{:});

      beam.angles = p.Results.angles;
    end
  end

  methods % Getters/setters
    function beam = set.angles(beam, val)
      assert(isnumeric(val) && numel(val) == 2 ...
          && min(val) >= 0.0 && max(val) <= pi, ...
          'angles must be 2 element numeric in range [0, pi]');
      beam.angles = sort(val(:).');
    end
  end
end
