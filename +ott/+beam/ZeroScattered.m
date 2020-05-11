classdef ZeroScattered < ott.beam.Scattered
% Scattered beam produced by zero-scattering methods.
% Inherits from :class:`ott.beam.Scattered`.
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
%   - incident_beam    -- The incident beam
%   - particle         -- The scattering method for force calculation
%   - type             -- Scattered type ('scattered', 'total' or 'internal')
%
% Dependent properties
%   - total_beam      -- Instance of the beam with total type
%   - scattered_beam  -- Instance of the beam with scattered type
%
% Scattered methods
%   - totalField      -- Calculate the total field beam
%   - scatteredField  -- Calculate the scattered field beam
%
% Methods deferred to `particle`
%   - force           -- Calls the particles force method
%   - torque          -- Calls the particles force method
%   - forcetorque     -- Calls the particles force method
%
% Methods deferred to `incident_beam`
%   - efieldInternal    -- Called by efield
%   - hfieldInternal    -- Called by hfield
%   - efarfieldInternal -- Called by efarfield
%   - hfarfieldInternal -- Called by hfarfield
%   - getBeamPower      -- get method called by dependent property power
%
% For details on methods, see :class:`Beam`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    particle           % The scattering method for force calculation
  end

  methods
    function beam = ZeroScattered(type, ibeam, particle, varargin)
      % Construct a new zero-scattered particle instance
      %
      % Usage
      %   beam = ZeroScattered(type, beam, particle, ...)
      %
      % Parameters
      %   - type (enum) -- Type of scattered beam.
      %   - beam (ott.beam.Beam) -- Incident beam.
      %   - particle (ott.scat.utils.Particle) -- Scattering method.
      %
      % Other parameters are passed to :class:`Scattered`.

      beam = beam@ott.beam.Scattered(type, 'incident_beam', ibeam, ...
          varargin{:});
      beam.particle = particle;
    end

    function force = force(beam, varargin)
      % Calculate force on a particle.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Uses the particle's force method for calculation.
      %
      % Usage
      %   force = beam.force(...)

      force = -beam.particle.force(beam.incident_beam, varargin{:});
    end

    function torque = torque(beam, varargin)
      % Calculate torque on a particle.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Uses the particle's torque method for calculation.
      %
      % Usage
      %   torque = beam.torque(...)

      torque = -beam.particle.torque(beam.incident_beam, varargin{:});
    end

    function [force, torque] = forcetorque(beam, varargin)
      % Calculate change in linear momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Uses the particle's forcetorque method for calculation.
      %
      % Usage
      %   [force, torque] = beam.forcetorque(...)

      [force, torque] = beam.particle.forcetorque(...
          beam.incident_beam, varargin{:});

      % Flip sign
      force = -force;
      torque = -torque;
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, xyz, varargin)
      % Get beam field (zeros or incident_beam)

      if strcmpi(beam.type, 'scattered')
        E = ott.utils.FieldVector(xyz, zeros(size(xyz)), 'cartesian');
      else
        % Defer to incident_beam
        E = beam.incident_beam.efieldInternal(xyz, varargin{:});
      end
    end

    function H = hfieldInternal(beam, xyz, varargin)
      % Get beam field (zeros or incident_beam)

      if strcmpi(beam.type, 'scattered')
        H = ott.utils.FieldVector(xyz, zeros(size(xyz)), 'cartesian');
      else
        % Defer to incident_beam
        H = beam.incident_beam.hfieldInternal(xyz, varargin{:});
      end
    end

    function E = efarfieldInternal(beam, rtp, varargin)
      % Get beam field (zeros or incident_beam)

      if strcmpi(beam.type, 'scattered')
        E = ott.utils.FieldVector(rtp, zeros(size(rtp)), 'spherical');
      else
        % Defer to incident_beam
        E = beam.incident_beam.efarfieldInternal(rtp, varargin{:});
      end
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Get beam field (zeros or incident_beam)

      if strcmpi(beam.type, 'scattered')
        H = ott.utils.FieldVector(rtp, zeros(size(rtp)), 'spherical');
      else
        % Defer to incident_beam
        H = beam.incident_beam.hfarfieldInternal(xyz, varargin{:});
      end
    end

    function p = getBeamPower(beam)
      % Get beam power (incident_beam or 0.0)

      if strcmpi(beam.type, 'scattered')
        p = 0.0;
      else
        % Defer to incident_beam
        p = beam.incident_beam.getBeamPower();
      end
    end
  end

  methods % Getters/setters
    function beam = set.particle(beam, val)
      assert(isa(val, 'ott.scat.utils.Particle'), ...
          'Particle should be a ott.scat.utils.Particle');
      assert(~isa(val, 'ott.scat.utils.BeamForce'), ...
          'Particle should not inherit from ott.scat.utils.BeamForce');
      beam.particle = val;
    end
  end
end
