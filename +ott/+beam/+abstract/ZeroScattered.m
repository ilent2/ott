classdef ZeroScattered < ott.beam.properties.ZeroScattered ...
    & ott.beam.abstract.ScatteredStem
% Describes a scattered beam with zero scattering.
% Inherits from :class:`ScatteredStem` and
% :class:`ott.beam.properties.ZeroScattered`.
%
% Zero-scattering methods do not scatter any of the incident beam, they
% are zero-order approximations.  The scattering methods support force
% calculation directly from the incident beam.  Calculating forces from
% the incident and scattered beam would produces no force.
%
% This class encapsulates the incident beam and the scattering method.
% Beam visualisation methods call the incident beam, while force
% calculation methods defer to the scattering method.
%
% Properties
%   - particle      -- Scattering method
%   - incident_beam -- Incident beam
%   - type          -- Default type used for casts/visualisation
%
% Dependent properties
%   - total_beam        -- The incident beam
%   - scattered_beam    -- An empty beam
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
    function beam = ZeroScattered(varargin)
      % Construct a new zero-scattered beam representation
      %
      % Usage
      %   beam = ZeroScattered(particle, incident_beam, ...)
      %
      % Parameters
      %   - particle (ott.scat.Particle) -- Scattering method
      %   - incident_beam (ott.beam.Beam) -- Incident beam
      %
      % Optional named arguments
      %   - type (enum) -- Default type for visualisations.
      %     Default: ``total``.

      beam = beam@ott.beam.properties.ZeroScattered(varargin{:});
    end

    function b = contains(beam, array_type)
      % Query if a array_type is contained in the array.
      %
      % Applies contains to the incident beam.
      %
      % Usage
      %   b = beam.contains(array_type)
      %
      % Parameters
      %   - array_type (enum) -- An array type, must be one of
      %     'array', 'coherent' or 'incoherent'.

      b = beam.incident_beam.contains(array_type);
    end

    function force = force(beam, varargin)
      % Calculate force on a particle.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Uses the particle's force method for calculation.
      %
      % Usage
      %   force = beam.force(...)

      p = inputParser;
      p.addOptional('other', [], ...
        @(x) isa(x, 'ott.beam.abstract.Beam') ...
        || isa(x, 'ott.scat.Particle'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if isempty(p.Results.other)
        force = -beam.particle.force(beam.incident_beam, unmatched{:});
      else
        force = beam.incident_beam.force(p.Results.other, unmatched{:});
      end
    end

    function torque = torque(beam, varargin)
      % Calculate torque on a particle.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Uses the particle's torque method for calculation.
      %
      % Usage
      %   torque = beam.torque(...)

      p = inputParser;
      p.addOptional('other', [], ...
        @(x) isa(x, 'ott.beam.abstract.Beam') ...
        || isa(x, 'ott.scat.Particle'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if isempty(p.Results.other)
        torque = -beam.particle.torque(beam.incident_beam, unmatched{:});
      else
        torque = beam.incident_beam.torque(p.Results.other, unmatched{:});
      end
    end

    function [force, torque] = forcetorque(beam, varargin)
      % Calculate change in linear momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Uses the particle's forcetorque method for calculation
      % unless a secondary particle/beam is provided.
      %
      % Usage
      %   [force, torque] = beam.forcetorque(...)
      %
      %   [force, torque] = beam.forcetorque(other, ...)

      p = inputParser;
      p.addOptional('other', [], ...
        @(x) isa(x, 'ott.beam.abstract.Beam') ...
        || isa(x, 'ott.scat.Particle'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if isempty(p.Results.other)
        [force, torque] = beam.particle.forcetorque(...
            beam.incident_beam, unmatched{:});

        % Flip sign
        force = -force;
        torque = -torque;
      else
        [force, torque] = beam.incident_beam.forcetorque(...
            p.Results.other, unmatched{:});
      end
    end
  end
end

