classdef HermiteGaussian < ott.beam.properties.Gaussian
% Properties of Hermite-Gaussian beams.
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist      -- Beam waist at focus
%   - mmode      -- Hermite mode order
%   - nmode      -- Hermite mode order
%   - power         -- Beam power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    mmode       % Hermite mode order
    nmode       % Hermite mode order
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.HermiteGaussian')
        args = ott.utils.addDefaultParameter('mmode', other.mmode, args);
        args = ott.utils.addDefaultParameter('nmode', other.nmode, args);
      end
      args = ott.beam.properties.Gaussian.likeProperties(other, args);
    end
  end
  
  methods
    function beam = HermiteGaussian(varargin)
      % Construct Hermite-Gaussian beam properties
      %
      % Usage
      %   beam = beam@ott.beam.properties.HermiteGaussian(...
      %       waist, mmode, nmode, ...)
      %   Parameters can also be passed as named arguments.
      
      p = inputParser;
      p.addOptional('waist', [], @isnumeric);
      p.addOptional('mmode', [], @isnumeric);
      p.addOptional('nmode', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      beam = beam@ott.beam.properties.Gaussian(p.Results.waist, unmatched{:});
      beam.mmode = p.Results.mmode;
      beam.nmode = p.Results.nmode;
    end
  end

  methods (Hidden)
    function b = isGaussian(beam)
      b = beam.mmode == 0 && beam.nmode == 0;
    end
  end

  methods % Getters/setters
    % mmode       % Hermite mode order
    % nmode       % Hermite mode order

    function beam = set.nmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'nmode must be numeric scalar');
      assert(round(val) == val, 'nmode must be integer');
      beam.nmode = val;
    end

    function beam = set.mmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'mmode must be numeric scalar');
      assert(round(val) == val, 'mmode must be integer');
      beam.mmode = val;
    end
  end
end
