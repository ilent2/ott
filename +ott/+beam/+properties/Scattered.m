classdef Scattered < ott.beam.properties.Beam
% Properties of scattered beams
% Inherits from :class:`Beam`.
%
% Properties
%   - beam_data       -- Beam object associated with the scattering
%   - incident_beam   -- The incident ray object (can be set to [])
%   - particle        -- Particle that caused the scattering (can be [])
%   - type            -- Type of beam ('scattered', 'total' or 'internal')
%
% Dependent properties
%   - omega           -- omega property of beam_data
%   - medium          -- medium property of beam_data
%   - power           -- Total beam power
%   - total_beam      -- Instance of the beam with total type
%   - scattered_beam  -- Instance of the beam with scattered type

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    beam_data           % Beam object associated with the scattering
    incident_beam       % The incident ray object (can be set to [])
    particle            %  Particle that caused the scattering (can be [])
    type                % Type of beam ('scattered', 'total' or 'internal')
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
        args = ott.utils.addDefaultParameter('type', other.type, args);
        args = ott.utils.addDefaultParameter(...
            'particle', other.particle, args);
        args = ott.utils.addDefaultParameter(...
            'incident_beam', other.incident_beam, args);
        args = ott.utils.addDefaultParameter(...
            'beam_data', other.beam_data, args);
      end
    end
  end

  methods
    function beam = Scattered(varargin)
      % Construct a new scattered beam
      %
      % Usage
      %   beam = beam@ott.beam.properties.Scattered(type, beam_data, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Optional parameters
      %   - particle ([]|ott.scat.Particle) -- Particle which caused
      %     the scattering.
      %
      %   - incident_beam ([]|ott.beam.Beam) -- Incident beam.

      p = inputParser;
      p.addOptional('type', [], ...
          @(x) any(strcmpi(x, {'incident', 'scattered', 'total'})));
      p.addOptional('beam_data', [], @(x) isa(x, 'ott.beam.Beam'));
      p.addParameter('particle', []);
      p.addParameter('incident_beam', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Beam(unmatched{:});
      beam = beam.setType(p.Results.type);
      beam.incident_beam = p.Results.incident_beam;
      beam.beam_data = p.Results.beam_data;
      beam.particle = p.Results.particle;
    end

    function beam = setType(beam, val)
      % Set the beam type paramter (without raising a warning)
      %
      % Usage
      %   beam = beam.setType(val);

      ott.utils.nargoutCheck(beam, nargout);

      S = warning('off', 'ott:beam:abstract:Scattered:type_change');
      beam.type = val;
      warning(S);
    end
  end

  methods % Getters/setters

    function beam = set.beam_data(beam, val)
      assert(isa(val, 'ott.beam.Beam'), ...
        'Beam data must be a Beam object');
      beam.beam_data = val;
    end

    function beam = set.incident_beam(beam, val)
      assert(isempty(val) || isa(val, 'ott.beam.Beam'), ...
        'Incident beam must be empty or a Beam object');
      beam.incident_beam = val;
    end

    function beam = set.particle(beam, val)
      assert(isempty(val) || isa(val, 'ott.scat.Particle'), ...
          'Particle must be empty of a scat.Particle object');
      beam.particle = val;
    end

    function beam = set.type(beam, val)

      % Check value
      assert(any(strcmpi(val, {'scattered', 'total', 'internal'})), ...
          'type must be ''scattered'', ''total'' or ''internal''');

      % Warn user they may be doing the wrong thing
      warning('ott:beam:abstract:Scattered:type_change', ...
        ['Changing the type property doesnt change the type', newline, ...
        'Consider using beam.total_beam or beam.scattered_beam instead.']);

      beam.type = val;
    end

    function omega = get.omega(beam)
      omega = beam.beam_data.omega;
    end
    function beam = set.omega(beam, val)

      % Don't change incident beam or scattering method
      warning('ott:beam:properties:Scattered:change_property', ...
          'Only changing property of scattered beam');

      beam.beam_data.omega = val;
    end

    function medium = get.medium(beam)
      medium = beam.beam_data.medium;
    end
    function beam = set.medium(beam, val)

      % Don't change incident beam or scattering method
      warning('ott:beam:properties:Scattered:change_property', ...
          'Only changing property of scattered beam');

      beam.beam_data.medium = val;
    end

    function power = get.power(beam)
      if strcmpi(beam.type, 'scattered')
        power = beam.total_beam.power;
      else
        power = beam.beam_data.power;
      end
    end
    function beam = set.power(beam, val)
      error('Changing scattered beam power not supported');
    end

    function tbeam = get.total_beam(beam)
      if strcmpi(beam.type, 'scattered')
        assert(~isempty(beam.incident_beam), ...
            'Need incident beam for conversion');
        tbeam = beam;
        tbeam.beam_data = beam.beam_data + beam.incident_beam;
        tbeam = tbeam.setType('total');
      elseif strcmpi(beam.type, 'total')
        tbeam = beam;
      else
        error('Unable to convert to specified type');
      end
    end

    function tbeam = get.scattered_beam(beam)
      if strcmpi(beam.type, 'scattered')
        tbeam = beam;
      elseif strcmpi(beam.type, 'total')
        assert(~isempty(beam.incident_beam), ...
            'Need incident beam for conversion');
        tbeam = beam;
        tbeam.beam_data = beam.beam_data - beam.incident_beam;
        tbeam = tbeam.setType('scattered');
      else
        error('Unable to convert to specified type');
      end
    end
  end
end
