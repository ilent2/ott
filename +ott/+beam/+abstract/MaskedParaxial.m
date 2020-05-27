classdef MaskedParaxial < ott.beam.abstract.MaskedFarfield ...
    & ott.beam.properties.FarfieldMapping
% Describes a composite beam with a masked paraxial field.
% Inherits from :class:`MaskedFarfield``.
%
% This beam is a container for another beam and applies a mask to
% to the paraxial far-field.  The beam can be cast to another type
% for calculating the near-field or far-field.
%
% Static methods
%   - TopHat        -- Creates a top-hat profile
%   - Annular       -- Creates a annular profile
%
% Properties
%   - beam          -- The masked beam (any ott.beam.Beam)
%   - mask          -- Function defining the mask
%   - mapping       -- Mapping for far-field to paraxial conversion
%   - direction     -- Mask direction (positive or negative hemisphere)
%
% Supported casts
%   - Beam
%   - Bsc
%   - abstract.InterpParaxial
%   - abstract.InterpFarfield -- (Inherited)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function beam = TopHat(radius, masked_beam, varargin)
      % Construct a top-hat beam
      %
      % Usage
      %   beam = MaskedParaxial.TopHat(radius, masked_beam, ...)

      mask = @(xy) vecnorm(xy) <= radius;
      beam = ott.beam.abstract.MaskedParaxial(mask, masked_beam, varargin{:});
    end

    function beam = Annular(radii, masked_beam, varargin)
      % Construct a paraxial masked annular beam
      %
      % Usage
      %   beam = MaskedParaxial.Annular(radii, masked_beam, ...)

      mask = @(xy) vecnorm(xy) >= radii(1) && vecnorm(xy) <= radii(2);
      beam = ott.beam.abstract.MaskedParaxial(mask, masked_beam, varargin{:});
    end
  end

  methods
    function beam = MaskedParaxial(mask, masked_beam, varargin)
      % Construct a masked paraxial beam
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
      %   - mapping (enum) -- 'sintheta' or 'tantheta'.  Default: 'tantheta'
      %   - hemisphere (enum) -- 'pos' or 'neg'.  Default: 'pos'

      p = inputParser;
      p.addParameter('mapping', 'tantheta');
      p.addParameter('hemisphere', 'pos');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.abstract.MaskedFarfield(...
          mask, masked_beam, unmatched{:});
      beam = beam@ott.beam.properties.FarfieldMapping(...
          'mapping', p.Results.mapping, 'hemisphere', p.Results.hemisphere);
    end
  end

  methods (Hidden)
    function E = efarfieldInternal(beam, rtp, varargin)
      % Calculate the internal beam far field and then mask

      % Calculate field at locations
      E = beam.beam.efarfield(rtp, varargin{:});

      % Compute mask
      xy = beam.farfield2paraxial(rtp, 'mapping', beam.mapping, ...
          'direction', beam.hemisphere);
      mask = beam.mask(xy);

      % Zero fields outside mask
      vrtp = E.vrtp(:, ~mask) = 0.0;
      E = ott.utils.FieldVector(E.locations, vrtp, 'spherical');
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Calculate the internal beam far field and then mask

      % Calculate field at locations
      H = beam.beam.efarfield(rtp, varargin{:});

      % Compute mask
      xy = beam.farfield2paraxial(rtp, 'mapping', beam.mapping, ...
          'direction', beam.hemisphere);
      mask = beam.mask(xy);

      % Zero fields outside mask
      vrtp = H.vrtp(:, ~mask) = 0.0;
      H = ott.utils.FieldVector(H.locations, vrtp, 'spherical');
    end
  end
end
