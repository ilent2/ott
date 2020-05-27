classdef MaskedFarfield < ott.beam.properties.Beam ...
    & ott.beam.abstract.CastNearfield
% Describes composite beam with a masked far-field.
% Inherits from :class:`CastNearfield` and :class:`ott.beam.properties.Beam`.
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
    beam
    mask
  end

  properties (Dependent)
    omega
    medium
    power
  end

  methods (Static)
    function beam = TopHat(angle, masked_beam, varargin)
      % Construct a top-hat beam
      %
      % Usage
      %   beam = MaskedParaxial.TopHat(angle, masked_beam, ...)

      mask = @(tp) tp(1, :) <= angle;
      beam = ott.beam.abstract.MaskedFarfield(mask, masked_beam, varargin{:});
    end

    function beam = Annular(angles, masked_beam, varargin)
      % Construct a paraxial masked annular beam
      %
      % Usage
      %   beam = MaskedParaxial.Annular(angles, masked_beam, ...)

      mask = @(tp) tp(1, :) >= angles(1) & tp(1, :) <= angles(2);
      beam = ott.beam.abstract.MaskedFarfield(mask, masked_beam, varargin{:});
    end
  end

  methods
    function beam = MaskedFarfield(mask, masked_beam, varargin)
      % Construct a masked far-field beam
      %
      % Usage
      %   beam = MaskedFarfield(mask, masked_beam, ...)
      %
      % Parameters
      %   - mask (function_handle) -- Masking function.  Should take a
      %     single argument for the paraxial xy position [x; y].
      %
      %   - masked_beam -- Beam to use for paraxial field calculations.

      beam = beam@ott.beam.properties.Beam(varargin{:});
      beam.mask = mask;
      beam.beam = masked_beam;
    end
  end

  methods (Hidden)
    function E = efarfieldInternal(beam, rtp, varargin)
      % Calculate the internal beam far field and then mask

      % Calculate field at locations
      E = beam.beam.efarfield(rtp, varargin{:});

      % Compute mask
      mask = beam.mask(rtp(2:3, :));

      % Zero fields outside mask
      vrtp = E.vrtp;
      vrtp(:, ~mask) = 0.0;
      E = ott.utils.FieldVector(E.locations, vrtp, 'spherical');
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Calculate the internal beam far field and then mask

      % Calculate field at locations
      H = beam.beam.efarfield(rtp, varargin{:});

      % Compute mask
      mask = beam.mask(rtp(2:3, :));

      % Zero fields outside mask
      vrtp = H.vrtp;
      vrtp(:, ~mask) = 0.0;
      H = ott.utils.FieldVector(H.locations, vrtp, 'spherical');
    end
  end

  methods % Getters/setters
    function p = get.power(beam)
      % How should we do this?  We did something fancy in the old TopHat?
      error('Not yet implemented');
    end

    function m = get.medium(beam)
      m = beam.beam.medium;
    end
    function m = get.omega(beam)
      m = beam.beam.omega;
    end
  end
end
