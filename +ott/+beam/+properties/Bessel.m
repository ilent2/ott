classdef Bessel < ott.beam.properties.Material
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

  methods (Static)
    function args = likeProperties(other, args)
      if isa(other, 'ott.beam.properties.Bessel')
        args = ott.utils.addDefaultParameter('angle', other.angle, args);
      end
      args = ott.beam.properties.Material.likeProperties(other, args);
    end
  end

  methods
    function beam = Bessel(varargin)
      % Construct Bessel properties
      %
      % Usage
      %   beam = beam@Bessel(angle, ...)
      %   Parameters can also be passed as named inputs.

      p = inputParser;
      p.addOptional('angle', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Material(unmatched{:});
      beam.angle = p.Results.angle;
    end
  end

  methods % Getters/setters
    function beam = set.angle(beam, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0.0 && val <= pi, ...
          'angle must be numeric scalar between 0 and pi');
      beam.angle = val;
    end
  end
end
