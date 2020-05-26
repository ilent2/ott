classdef Scattered < ott.beam.Scattered ...
    & ott.beam.abstract.Beam
% Represents a scattered beam.
% Inherits from :class:`Beam` and :class:`ott.beam.Scattered`.
%
% The scattered beam class is a container which encapsulates a scattered
% beams type and the incident beam that cause the scattering.  Scattered
% beams behave like other beams, most properties are dependent on the
% internal beam_data property.
%
% Supported casts
%   - Beam            -- Default cast to Beam, uses Scattered
%   - Scattered       -- Cast the Scattered object but not the beam_data

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function beam = Scattered(type, beam, varargin)
      % Construct a new scattered beam representation
      %
      % Usage
      %   beam = Scattered(type, beam, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - type (enum) -- Type of scattered beam.
      %     Either 'scattered', 'total' or 'internal'.
      %
      %   - beam (ott.beam.Beam) -- The beam data.
      %
      % Optional named parameters
      %   - incident_beam ([]|Beam) -- Incident beam or emtpy.
      %     Default: ``[]``.
      %
      %   - particle ([]|ott.scat.utils.Particle) -- The particle that
      %     caused the scattering or empty.  Default ``[]``.

      beam = beam@ott.beam.properties.Scattered(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to a Scattered beam
      beam = ott.beam.Scattered(varargin{:});
    end

    function beam = ott.beam.Scattered(beam, varargin)
      % Cast the beam to a scattered beam, leaving internal data untouched.
      beam = castHelper(@ott.beam.Scattered.like, ...
          beam, varargin{:});
    end
  end
end
