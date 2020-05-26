classdef LaguerreGaussian < ott.beam.properties.Gaussian
% Properties of paraxial Laguerre-Gaussian beams.
% Inherits from :class:`ott.beam.properties.Gaussian`.
%
% Properties
%   - waist         -- Beam waist radius
%   - lmode         -- Azimuthal Laguerre mode order
%   - pmode         -- Radial Laguerre mode order
%   - power         -- Beam power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    lmode       % Azimuthal Laguerre mode order
    pmode       % Radial Laguerre mode order
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.LaguerreGaussian')
        args = ott.utils.addDefaultParameter('lmode', other.lmode, args);
        args = ott.utils.addDefaultParameter('pmode', other.pmode, args);
      end
      args = ott.beam.properties.Gaussian.likeProperties(other, args);
    end
  end

  methods
    function beam = LaguerreGaussian(varargin)
      % Construct LG beam properties

      p = inputParser;
      p.addOptional('waist', [], @isnumeric);
      p.addOptional('lmode', [], @isnumeric);
      p.addOptional('pmode', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Gaussian(p.Results.waist, unmatched{:});
      beam.lmode = p.Results.lmode;
      beam.pmode = p.Results.pmode;
    end
  end

  methods % Getters/setters
    % waist       % Beam waist at focus
    % lmode       % Azimuthal Laguerre mode order
    % pmode       % Radial Laguerre mode order

    function beam = set.lmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'lmode must be numeric scalar');
      assert(round(val) == val, ...
          'lmode must be integer');
      beam.lmode = val;
    end

    function beam = set.pmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'pmode must be numeric scalar');
      assert(round(val) == val && val >= 0, ...
          'pmode must be positive integer');
      beam.pmode = val;
    end
  end
end
