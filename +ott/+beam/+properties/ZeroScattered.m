classdef ZeroScattered < ott.beam.properties.ScatteredInterface
% Properties for zero-scattered beams.
% Inherits from :class:`ScatteredInterface`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    incident_beam  % Incident beam
    particle       % Scattering particle
  end

  properties (Dependent)
    omega               % omega property of beam_data
    medium              % medium property of beam_data
    power               % Total beam power
    total_beam          % Instance of the beam with total type
    scattered_beam      % Instance of the beam with scattered type
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Adds similar properties to an array
      if isa(other, 'ott.beam.properties.Scattered')
        args = ott.utils.addDefaultParameter(...
            'particle', other.particle, args);
        args = ott.utils.addDefaultParameter(...
            'incident_beam', other.incident_beam, args);
      end
      args = ott.beam.properties.ScatteredInterface.likeProperties(...
          other, args);
    end
  end

  methods
    function beam = ZeroScattered(varargin)
      % Construct zero-scattered properties
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

      p = inputParser;
      p.addOptional('particle', [], @(x) isa(x, 'ott.scat.Particle'));
      p.addOptional('incident_beam', [], @(x) isa(x, 'ott.beam.Beam'));
      p.addParameter('type', 'total');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.ScatteredInterface(...
          'type', p.Results.type, unmatched{:});
      beam.incident_beam = p.Results.incident_beam;
      beam.particle = p.Results.particle;
    end
  end

  methods % Getters/setters
    function beam = set.incident_beam(beam, val)
      assert(isa(val, 'ott.beam.Beam'), ...
        'Incident beam must be a Beam object');
      beam.incident_beam = val;
    end

    function beam = set.particle(beam, val)
      assert(isa(val, 'ott.scat.utils.ZeroScattered'), ...
          'particle must be a scat.utils.ZeroScattered object');
      beam.particle = val;
    end

    function sbeam = get.scattered_beam(beam)
      sbeam = ott.beam.abstract.Empty();
    end

    function tbeam = get.total_beam(beam)
      tbeam = beam.incident_beam;
    end

    function p = get.power(beam)
      p = beam.incident_beam.power;
    end
    function beam = set.power(beam, val)
      beam.incident_beam.power = val;
    end

    function p = get.omega(beam)
      p = beam.incident_beam.omega;
    end
    function beam = set.omega(beam, val)
      beam.incident_beam.omega = val;
    end

    function p = get.medium(beam)
      p = beam.incident_beam.medium;
    end
    function beam = set.medium(beam, val)
      beam.incident_beam.medium = val;
    end
  end
end

