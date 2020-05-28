classdef InceGaussian < ott.beam.properties.Gaussian
% Properties of Ince-Gaussian beams.
% Inherits from :class:`Gaussian`.
%
% Properties
%   - waist       -- Beam waist radius at focus
%   - lmode       -- Azimuthal mode number
%   - porder      -- Paraxial mode order.
%   - ellipticity -- Ellipticity of coordinates
%   - parity      -- Parity of beam ('even' or 'odd')
%   - power       -- Beam power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    lmode        % Azimuthal mode number
    porder       % Paraxial mode order.
    ellipticity  % Ellipticity of coordinates
    parity       % Parity of beam ('even' or 'odd')
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.InceGaussian')
        args = ott.utils.addDefaultParameter('lmode', other.lmode, args);
        args = ott.utils.addDefaultParameter('porder', other.porder, args);
        args = ott.utils.addDefaultParameter(...
            'ellipticity', other.ellipticity, args);
        args = ott.utils.addDefaultParameter('parity', other.parity, args);
      end
      args = ott.beam.properties.Gaussian.likeProperties(other, args);
    end
  end

  methods
    function beam = InceGaussian(varargin)
      % TODO: Docs

      p = inputParser;
      p.addOptional('waist', [], @(x) isnumeric(x) & ~isempty(x));
      p.addOptional('lmode', [], @(x) isnumeric(x) & ~isempty(x));
      p.addOptional('porder', [], @(x) isnumeric(x) & ~isempty(x));
      p.addOptional('parity', [], @(x) any(strcmpi(x, {'even', 'odd'})));
      p.addOptional('ellipticity', 1, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Gaussian(p.Results.waist, unmatched{:});
      beam.parity = p.Results.parity;
      beam.ellipticity = p.Results.ellipticity;
      beam.lmode = p.Results.lmode;
      beam.porder = p.Results.porder;
    end
  end

  methods (Hidden)
    function b = isGaussian(beam)
      % TODO: Implement
      error('Not yet implemented');
    end
  end

  methods % Getters/setters
    % lmode        % Azimuthal mode number
    % pmode        % Paraxial mode order.
    % ellipticity  % Ellipticity of coordinates
    % parity       % Parity of beam ('even' or 'odd')

    function beam = set.lmode(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'lmode must be numeric scalar');
      assert(round(val) == val, 'lmode must be integer');
      beam.lmode = val;
    end

    function beam = set.porder(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'porder must be numeric scalar');
      assert(round(val) == val, 'pmode must be integer');
      beam.porder = val;
    end

    function beam = set.ellipticity(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'ellipticity must be numeric scalar');
      beam.ellipticity = val;
    end

    function beam = set.parity(beam, val)
      assert(any(strcmpi(val, {'even', 'odd'})), ...
          'parity must be ''even'' or ''odd''');
      beam.parity = val;
    end
  end
end
