classdef Ray < ott.beam.abstract.PlaneWaveStem
% Specialisation describing abstract geometric ray beams.
% Inherits from :class:`ott.beam.properties.Ray` and :class:`Beam`.
%
% This class is only provided for consistency.  To create a Ray beam,
% it is recommended to create either a :class:`ott.beam.Ray` directly
% or create another beam and cast to a `ott.beam.Ray`, for example
%
% .. code:: matlab
%   beam = ott.beam.abstract.Gaussian(1.0);
%   rays = ott.beam.Ray(beam);
%
% Static methods
%   - like            -- Create a beam like another
%   - FromDirection   -- Construct from direction/polarisation
%   - DirectionSet    -- Construct a directionSet
%
% Dependent properties
%   - power      -- Power of ray
%
% Casts
%   - Beam       -- Uses Ray
%
% All other casts inherited from base.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    power
  end

  methods (Static)
    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Ray.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.Ray.likeProperties(other, varargin);
      beam = ott.beam.abstract.Ray(args{:});
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

      % Construct direction set
      directionSet = ott.beam.properties.PlaneWave.DirectionSet(...
          p.Reults.direction, p.Results.polarisation1);

      % Construct beam
      beam = ott.beam.abstract.Ray(...
          'origin', p.Results.origin, ...
          'directionSet', directionSet, ...
          'field', p.Results.field, ...
          unmatched{:});
    end
  end

  methods
    function varargout = visualise(beam, varargin)
      % Cast abstract.Ray to Ray and visualise
      %
      % Usage
      %   [...] = beam.visualise(...)

      beam = ott.beam.Ray(beam);
      beam.visualise(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to Ray
      beam = ott.beam.Ray(varargin{:});
    end
  end

  methods % Getters/setters
    function p = get.power(beam)
      p = beam.intensity;
    end
  end
end
