classdef MaskedNearfield3d < ott.beam.properties.MaskedBeam ...
    & ott.beam.abstract.CastFarfield
% Describes composite beam with a masked near-field.
% Inherits from :class:`CastFarfield` and :class:`ott.beam.properties.Beam`.
%
% Static methods
%   - TopHat        -- Creates a top-hat profile
%   - Annular       -- Creates a annular profile
%
% Properties
%   - beam          -- The masked beam (any ott.beam.Beam)
%   - mask          -- Function defining the mask
%
% All casts inherited from base.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    maxRadius
  end

  methods
    function beam = MaskedNearfield3d(varargin)
      % Construct a masked far-field beam
      %
      % Usage
      %   beam = MaskedNearfield3d(mask, masked_beam, ...)
      %
      % Parameters
      %   - mask (function_handle) -- Masking function.  Should take a
      %     single argument for the paraxial xy position [x; y].
      %
      %   - masked_beam -- Beam to use for paraxial field calculations.
      %
      %   - maxRadius (numeric) -- Maximum radius of mask.
      %     Default: ``Inf``.

      p = inputParser;
      p.addOptional('mask', [], @(x) isa(x, 'function_handle'));
      p.addOptional('masked_beam', [], @(x) isa(x, 'ott.beam.Beam'));
      p.addParameter('maxRadius', Inf);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.MaskedBeam(...
        'mask', p.Results.mask, ...
        'masked_beam', p.Results.masked_beam, ...
        unmatched{:});
      beam.maxRadius = p.Results.maxRadius;
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, xyz, varargin)
      % Calculate the internal beam near field and then mask

      % TODO: Support projection
      assert(all(vecnorm(xyz) <= beam.maxRadius), ...
          'all points must be maxRadius (for now)');

      % Calculate field at locations
      E = beam.masked_beam.efield(xyz, varargin{:});

      % Compute mask
      mask = beam.mask(xyz);

      % Zero fields outside mask
      vxyz = E.vxyz;
      vxyz(:, ~mask) = 0.0;
      E = ott.utils.FieldVector(xyz, vxyz, 'cartesian');
    end

    function H = hfieldInternal(beam, xyz, varargin)
      % Calculate the internal beam near field and then mask

      % TODO: Support projection
      assert(all(vecnorm(xyz) <= beam.maxRadius), ...
          'all points must be maxRadius (for now)');

      % Calculate field at locations
      H = beam.masked_beam.hfield(rtp, varargin{:});

      % Compute mask
      mask = beam.mask(xyz);

      % Zero fields outside mask
      vxyz = H.vxyz;
      vxyz(:, ~mask) = 0.0;
      H = ott.utils.FieldVector(xyz, vxyz, 'cartesian');
    end
  end
end
