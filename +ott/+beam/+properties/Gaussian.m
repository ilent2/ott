classdef Gaussian < ott.beam.properties.Beam
% Properties of a paraxial Gaussian beam.
% Inherits from :class:`ott.beam.properties.Beam`.
%
% Properties
%   - waist         -- Beam waist radius
%   - power         -- Beam power
%   - polarisation  -- Polarisation of the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    waist          % Beam waist radius
    polarisation   % Polarisation of the beam
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      args = ott.utils.addDefaultParameter('waist', other.waist, args);
      args = ott.utils.addDefaultParameter(...
          'polarisation', other.polarisation, args);
      args = ott.beam.properties.Beam.likeProperties(other, args);
    end
  end

  methods
    function beam = Gaussian(varargin)
      % Construct a new Gaussian beam
      %
      % Usage
      %   beam = beam@ott.beam.properties.Gaussian(waist, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist.
      %
      % Optional named arguments
      %   - polarisation (2 numeric) -- Polarisation of the beam.
      %     Default: ``[1, 1i]``.

      p = inputParser;
      p.addOptional('waist', [], @isnumeric);
      p.addParameter('polarisation', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Beam(unmatched{:});
      beam.polarisation = p.Results.polarisation;
      beam.waist = p.Results.waist;
    end
  end

  methods % Getters/setters
    function beam = set.waist(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'waist must be numeric scalar');
      beam.waist = val;
    end

    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && numel(val) == 2 && isvector(val), ...
          'polarisation must be 2 element numeric vector');
      beam.polarisation = val(:).';
    end
  end
end
