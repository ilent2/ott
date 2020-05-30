classdef MaskedParaxial < ott.beam.properties.MaskedBeam ...
    & ott.beam.properties.FarfieldMapping
% Properties of masked paraxial beams

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      args = ott.beam.properties.FarfieldMapping.likeProperties(other, args);
      args = ott.beam.properties.MaskedBeam.likeProperties(other, args);
    end
  end

  methods
    function beam = MaskedParaxial(varargin)
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
      %
      % Optional named arguments
      %   - mapping (enum) -- 'sin', 'tan' or 'theta'.  Default: 'tan'
      %   - hemisphere (enum) -- 'pos' or 'neg'.  Default: 'pos'

      p = inputParser;
      p.addOptional('mask', [], @(x) isa(x, 'function_handle'));
      p.addOptional('masked_beam', [], @(x) isa(x, 'ott.beam.Beam'));
      p.addParameter('mapping', 'tan');
      p.addParameter('hemisphere', 'pos');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.MaskedBeam(...
          'mask', p.Results.mask, ...
          'masked_beam', p.Results.masked_beam, ...
          unmatched{:});
      beam = beam@ott.beam.properties.FarfieldMapping(...
          'mapping', p.Results.mapping, 'hemisphere', p.Results.hemisphere);
    end
  end
end
