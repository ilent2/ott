classdef MaskedBeam < ott.beam.properties.Beam
% Properties of masked beams
% Inherits from :class:`Beam`.
%
% Properties
%   - masked_beam
%   - mask
%
% Dependent properties
%    - omega      -- omega property of masked_beam
%    - medium     -- medium property of masked_beam
%    - power      -- not yet implemented

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    masked_beam
    mask
  end

  properties (Dependent)
    omega
    medium
    power
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.MaskedBeam')
        args = ott.utils.addDefaultParameter('mask', other.mask, args);
        args = ott.utils.addDefaultParameter(...
            'masked_beam', other.masked_beam, args);
      end
      args = ott.beam.properties.Beam.likeProperties(other, args);
    end
  end


  methods
    function beam = MaskedBeam(varargin)
      % Construct masked paraxial properties
      %
      % Usage
      %   beam = MaskedParaxial(mask, masked_beam, ...)
      %
      % Parameters
      %   - mask (function_handle) -- Masking function.  Should take a
      %     single argument for the paraxial xy position [x; y].
      %
      %   - masked_beam -- Beam to use for paraxial field calculations.

      p = inputParser;
      p.addOptional('mask', [], @(x) isa(x, 'function_handle'));
      p.addOptional('masked_beam', [], @(x) isa(x, 'ott.beam.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Beam(unmatched{:});
      beam.mask = p.Results.mask;
      beam.masked_beam = p.Results.masked_beam;
    end
  end

  methods % Getters/setters
    function beam = set.mask(beam, val)
      assert(isa(val, 'function_handle'), 'mask must be function handle');
      beam.mask = val;
    end

    function beam = set.masked_beam(beam, val)
      assert(isa(val, 'ott.beam.Beam'), 'beam must be ott.beam.Beam');
      beam.masked_beam = val;
    end

    function p = get.power(beam)
      % How should we do this?  We did something fancy in the old TopHat?
      error('Not yet implemented');
    end

    function m = get.medium(beam)
      m = beam.masked_beam.medium;
    end
    function beam = set.medium(beam, val)
      beam.masked_beam.medium = val;
    end

    function m = get.omega(beam)
      m = beam.masked_beam.omega;
    end
    function beam = set.omega(beam, val)
      beam.masked_beam.omega = val;
    end
  end
end
