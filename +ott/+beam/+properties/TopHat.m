classdef TopHat < ott.beam.properties.Material
% Properties of a collimated Top-Hat beam.
%
% Properties
%   - radius        -- Radius of the top-hat
%   - polarisation  -- Polarisation of the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    radius          % Radius of the beam
    polarisation    % Polarisation of the beam
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.TopHat')
        args = ott.utils.addDefaultParameter('radius', other.radius, args);
        args = ott.utils.addDefaultParameter('polarisation', ...
            other.polarisation, args);
      end
      args = ott.beam.properties.Material.likeProperties(other, args);
    end
  end

  methods
    function beam = TopHat(varargin)
      % Construct TopHat properties
      %
      % Usage
      %   beam = beam@TopHat(radius, polarisation)
      %   Parameters can also be named.

      p = inputParser;
      p.addOptional('radius', [], @isnumeric);
      p.addOptional('polarisation', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Material(unmatched{:});
      beam.radius = p.Results.radius;
      beam.polarisation = p.Results.polarisation;
    end
  end

  methods % Getters/setters
    function beam = set.radius(beam, val)
      assert(isnumeric(val) && isscalar(val) && val > 0.0, ...
          'radius must be positive numeric scalar');
      beam.radius = val;
    end
    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'polarisation must be 2 element numeric vector');
      beam.polarisation = val;
    end
  end
end
