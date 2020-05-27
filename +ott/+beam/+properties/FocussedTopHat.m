classdef FocussedTopHat
% Properties of a focussed Top-hat beam
%
% Properties
%   - angle           -- Angle of top-hat edge

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    angle
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.FocussedTopHat')
        args = ott.utils.addDefaultParameter('angle', other.angle, args);
      end
    end
  end

  methods
    function beam = FocussedTopHat(varargin)
      % Construct Focussed-TopHat properties
      %
      % Usage
      %   beam = beam@FocussedTopHat(angle)

      p = inputParser;
      p.addOptional('angle', [], @isnumeric);
      p.parse(varargin{:});

      beam.angle = p.Results.angle;
    end
  end

  methods % Getters/setters
    function beam = set.angle(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'angle must be numeric scalar');
    end
  end
end
