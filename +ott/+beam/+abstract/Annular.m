classdef Annular < ott.beam.properties.Annular ...
    & ott.beam.properties.Material ...
    & ott.beam.abstract.CastBoth ...
    & ott.beam.properties.VariablePower
% Abstract description of annular beams.
% Inherits from :class:`ott.beam.properties.Annular` and
% :class:`ott.beam.properties.Material` and :class:`CastBoth`.
%
% Properties
%   - angles
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
      args = ott.beam.properties.Annular.likeProperties(other, args);
      args = ott.beam.properties.VariablePower.likeProperties(other, args);
      args = ott.beam.abstract.CastBoth.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Annular.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.Annular.likeProperties(other, varargin);
      beam = ott.beam.abstract.Annular(args{:});
    end
  end

  methods
    function beam = Annular(varargin)
      % Constructs a new representation of a annular beam
      %
      % Usage
      %   beam = Annular(angles, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Optional named arguments
      %   - power (numeric) -- beam power.  Default: ``1.0``.

      p = inputParser;
      p.addOptional('angles', [], @isnumeric);
      p.addParameter('power', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Annular(p.Results.angles);
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

      p = inputParser;
      p.addParameter('masked_beam', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % TODO: Multiple beams and type check
      % TODO: Use beam power (needs to be scaled by mask area)

      % Get beam to mask
      masked_beam = p.Results.masked_beam;
      if isempty(masked_beam)
        masked_beam = ott.beam.abstract.UniformFarfield([1;0]);
      end

      beam = ott.beam.abstract.MaskedFarfield.Annular(beam.angles, ...
          masked_beam, unmatched{:});
    end
  end
end
