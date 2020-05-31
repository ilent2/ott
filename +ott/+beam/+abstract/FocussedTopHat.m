classdef FocussedTopHat < ott.beam.properties.FocussedTopHat ...
    & ott.beam.properties.Material ...
    & ott.beam.abstract.CastBoth ...
    & ott.beam.properties.VariablePower
% Abstract description of a focussed top-hat beam.
% Inherits from :class:`ott.beam.properties.FocussedTopHat` and
% :class:`ott.beam.properties.Material` and :class:`CastBoth`.
%
% Properties
%   - angle
%   - power
%
% Casts
%   - Beam           -- Overloaded, Uses abstract.MaskedFarfield
%   - abstract.MaskedFarfield
%
% Additional casts inherited from base.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.properties.FocussedTopHat.likeProperties(other, args);
      args = ott.beam.properties.Material.likeProperties(other, args);
      args = ott.beam.properties.VariablePower.likeProperties(other, args);
      args = ott.beam.abstract.CastBoth.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = FocussedTopHat.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.FocussedTopHat.likeProperties(other, varargin);
      beam = ott.beam.abstract.FocussedTopHat(args{:});
    end
  end

  methods
    function beam = FocussedTopHat(varargin)
      % Constructs a new representation of a focussed top hat beam
      %
      % Usage
      %   beam = FocussedTopHat(angle, ...)
      %
      % Optional named arguments
      %   - power (numeric) -- beam power.  Default: ``1.0``.

      p = inputParser;
      p.addOptional('angle', [], @isnumeric);
      p.addParameter('power', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.FocussedTopHat(p.Results.angle);
      beam = beam@ott.beam.properties.Material(unmatched{:}, ...
          'power', p.Results.power);
    end

    function beam = ott.beam.Beam(varargin)
      % Cast the beam to a MaskedFarfield
      beam = ott.beam.abstract.MaskedFarfield(varargin{:});
    end

    function beam = ott.beam.abstract.MaskedFarfield(beam, varargin)
      % Construct a masked far-field beam
      %
      % Usage
      %   beam = MaskedFarfield(beam, masked_beam, ...)
      %   For additional arguments see :class:`MaskedFarfield`.
      %
      % Optional arguments
      %   - masked_beam -- Beam to use for mask.
      %     Default: ``UniformFarfield([1;0])``.

      assert(isa(beam, 'ott.beam.abstract.FocussedTopHat'), ...
          'First argument must be a abstract.FocussedTopHat');
      ott.utils.nargoutCheck(beam, nargout);

      p = inputParser;
      p.addParameter('masked_beam', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get beam to mask
      masked_beam = p.Results.masked_beam;
      if isempty(masked_beam)
        masked_beam = ott.beam.abstract.UniformFarfield([1;0]);
        masked_beam.power = 4*pi./(2*pi*beam.angle) .* beam.power;
      end

      beam = ott.beam.abstract.MaskedFarfield.TopHat(beam.angle, ...
          masked_beam, unmatched{:});
    end
  end
end
