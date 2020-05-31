classdef PlaneWave < ott.beam.abstract.PlaneWaveStem ...
    & ott.beam.properties.InfinitePower
% Abstract representation of a plane wave beam.
% Inherits from :class:`Beam`
% and :class:`ott.beam.properties.PlaneWaveScalar`.
%
% This class creates a single abstract Plane Wave beam.  For more
% efficient arrays of plane wave beams, use :class:`ott.beam.PlaneWave`
% directly.
%
% Static methods
%   - like            -- Create a beam like another
%   - FromDirection   -- Construct from direction/polarisation
%   - DirectionSet    -- Construct a directionSet
%
% Casts
%   - Beam      -- (Sealed) Uses PlaneWave
%
% All other casts inherited from base.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = PlaneWave.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.PlaneWave.likeProperties(other, varargin);
      beam = ott.beam.abstract.PlaneWave(args{:});
    end

    function beam = FromDirection(varargin)
      % Construct beam from direction/polarisation vectors.
      %
      % Usage
      %   beam = FromDirection(origin, direction, polarisation, field)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - origin (3xN numeric) -- Origin (for phase offset) of wave.
      %   - direction (3xN numeric) -- Propagation direction of wave.
      %   - polarisation (3xN numeric) -- Primary polarisation direction.
      %   - field (2xN numeric) -- Field in two polarisation directions.
      %
      % Additional parameters passed to base.

      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('direction', [], @isnumeric);
      p.addOptional('polarisation1', [], @isnumeric);
      p.addOptional('field', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct direction set
      directionSet = ott.beam.properties.PlaneWave.DirectionSet(...
          p.Results.direction, p.Results.polarisation1);

      % Construct beam
      beam = ott.beam.abstract.PlaneWave(...
          'origin', p.Results.origin, ...
          'directionSet', directionSet, ...
          'field', p.Results.field, ...
          unmatched{:});
    end
  end

  methods (Sealed)
    function beam = ott.beam.Beam(varargin)
      % Cast to PlaneWave
      beam = ott.beam.PlaneWave(varargin{:});
    end
  end
end
