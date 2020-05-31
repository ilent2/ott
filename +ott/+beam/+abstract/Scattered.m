classdef Scattered < ott.beam.abstract.ScatteredStem ...
    & ott.beam.properties.Scattered
% Represents a scattered beam.
% Inherits from :class:`CastBoth` and :class:`ott.beam.properties.Scattered`.
%
% This class encapsulate two beam objects: the incident beam and the
% scattered beam.  Scattered beams can describe either the total field
% or only the scattered field.  This class provides an uniform encapsulation
% of beams in either of these types and methods to switch between the
% two when the incident field is known.
%
% The incident beam assigned to this object need not be the same incident
% beam involved in the scattering calculation.  This can facilitate more
% efficient scattering calculations.
%
% This class only represents external fields.
% To change the behaviour of the total_beam and scattered_beam properties,
% overload the `plus` and `minus` operators for the internal beams.
%
% Properties
%   - incident_beam -- Incident beam
%   - data          -- Scattered beam data
%   - data_type     -- Scattered beam type
%   - type          -- Default type used for casts/visualisations
%
% Dependent properties
%   - total_beam      -- Total-field beam (incident+scattered)
%   - scattered_beam  -- Scattered-field beam
%
% Methods
%   - force         -- Simplified force calculation
%   - torque        -- Simplified torque calculation
%   - forcetorque   -- Simplified force/torque calculation
%
% Casts
%   - Beam          -- Retrieves the total_beam or scattered_beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function beam = Scattered(varargin)
      % Construct a new scattered beam representation
      %
      % Usage
      %   beam = Scattered(data_type, data, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - data_type (enum) -- Type of scattered beam.
      %     Either 'scattered' or 'total'.
      %
      %   - data (ott.beam.Beam) -- The beam data.
      %
      % Optional named parameters
      %   - incident_beam ([]|Beam) -- Incident beam or empty.
      %     Default: ``[]``.
      %
      %   - type (enum) -- Default type used for casts and visualisations.
      %     Default: ``data_type``.

      beam = beam@ott.beam.properties.Scattered(varargin{:});
    end

    function b = contains(beam, array_type)
      % Query if a array_type is contained in the array.
      %
      % Casts the beam and applies the contains method.
      %
      % Usage
      %   b = beam.contains(array_type)
      %
      % Parameters
      %   - array_type (enum) -- An array type, must be one of
      %     'array', 'coherent' or 'incoherent'.

      b = ott.beam.Beam(beam).contains(array_type);
    end

    function varargout = force(ibeam, varargin)
      % Calculate change in linear momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   force = ibeam.force(other, ...)
      %
      % If other is not supplied, uses the incident beam for the calculation.

      p = inputParser;
      p.addOptional('other', [], ...
        @(x) isa(x, 'ott.beam.abstract.Beam') ...
        || isa(x, 'ott.scat.Particle'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      other = p.Results.other;
      ibeam = ibeam.total_beam;

      if isempty(other)
        if isempty(ibeam.incident_beam)
          error('Must supply other or have valid incident_beam');
        end

        other = ibeam.incident_beam;
      end

      [varargout{1:nargout}] = ibeam.force(other, unmatched{:});
    end

    function varargout = torque(ibeam, varargin)
      % Calculate change in angular momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   torque = ibeam.torque(other, ...)
      %
      % If other is not supplied, uses the incident beam for the calculation.

      p = inputParser;
      p.addOptional('other', [], ...
        @(x) isa(x, 'ott.beam.abstract.Beam') ...
        || isa(x, 'ott.scat.Particle'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      other = p.Results.other;
      ibeam = ibeam.total_beam;

      if isempty(other)
        if isempty(ibeam.incident_beam)
          error('Must supply other or have valid incident_beam');
        end

        other = ibeam.incident_beam;
      end

      [varargout{1:nargout}] = ibeam.torque(other, unmatched{:});
    end

    function varargout = forcetorque(ibeam, varargin)
      % Calculate change in momentum between beams.
      %
      % Usage
      %   [force, torque] = ibeam.forcetorque(other, ...)
      %   Returns 3xN matrices for the force and torque in Cartesian
      %   coordinates.
      %
      % If other is not supplied, uses the incident beam for the calculation.
      %
      % Parameters
      %   - other (Beam|scat.Scatter) -- A beam to compare the force
      %     with or a particle with a scatter method.
      %
      %   - position (3xN numeric) -- Distance to translate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Default: ``[]``.
      %
      %   - rotation (3x3N numeric) -- Angle to rotate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Inverse rotation is applied to scattered beam, effectively
      %     rotating the particle.
      %     Default: ``[]``.

      p = inputParser;
      p.addOptional('other', [], ...
        @(x) isa(x, 'ott.beam.abstract.Beam') ...
        || isa(x, 'ott.scat.Particle'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      other = p.Results.other;
      ibeam = ibeam.total_beam;

      if isempty(other)
        if isempty(ibeam.incident_beam)
          error('Must supply other or have valid incident_beam');
        end

        other = ibeam.incident_beam;
      end

      [varargout{1:nargout}] = ibeam.forcetorque(other, unmatched{:});
    end
  end
end
