classdef Scattered < ott.beam.properties.ScatteredInterface
% Properties of scattered beams
% Inherits from :class:`ScatteredInterface`.
%
% Properties
%   - incident_beam -- Incident beam
%   - data          -- Scattered beam data
%   - data_type     -- Scattered beam type
%   - type          -- Default type used for casts/visualisations
%
% Dependent properties
%   - omega           -- omega property of beam_data
%   - medium          -- medium property of beam_data
%   - power           -- Total beam power
%   - total_beam      -- Total-field beam (incident+scattered)
%   - scattered_beam  -- Scattered-field beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    incident_beam  % Incident beam
    data           % Scattered beam data
    data_type      % Scattered beam type
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
            'data_type', other.data_type, args);
        args = ott.utils.addDefaultParameter(...
            'data', other.data, args);
        args = ott.utils.addDefaultParameter(...
            'incident_beam', other.incident_beam, args);
      end
      args = ott.beam.properties.ScatteredInterface.likeProperties(...
          other, args);
    end
  end

  methods
    function beam = Scattered(varargin)
      % Construct a new scattered beam
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

      p = inputParser;
      p.addOptional('data_type', [], ...
          @(x) any(strcmpi(x, {'scattered', 'total'})));
      p.addOptional('data', [], @(x) isa(x, 'ott.beam.Beam'));
      p.addParameter('type', []);
      p.addParameter('incident_beam', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if ~isempty(p.Results.type)
        type = p.Results.type;
      else
        type = p.Results.data_type;
      end

      beam = beam@ott.beam.properties.ScatteredInterface(...
          'type', type, unmatched{:});

      beam.data_type = p.Results.data_type;
      beam.data = p.Results.data;
      beam.incident_beam = p.Results.incident_beam;
    end
  end

  methods % Getters/setters

    function beam = set.data(beam, val)
      assert(isa(val, 'ott.beam.Beam'), ...
        'Beam data must be a Beam object');
      beam.data = val;
    end

    function beam = set.incident_beam(beam, val)
      assert(isempty(val) || isa(val, 'ott.beam.Beam'), ...
        'Incident beam must be empty or a Beam object');
      beam.incident_beam = val;
    end

    function beam = set.data_type(beam, val)
      assert(any(strcmpi(val, {'scattered', 'total'})), ...
          'type must be ''scattered'' or ''total''');
      beam.data_type = val;
    end

    function omega = get.omega(beam)
      omega = beam.data.omega;
    end
    function beam = set.omega(beam, val)
      beam.data.omega = val;
      beam.incident_beam.omega = val;
    end

    function medium = get.medium(beam)
      medium = beam.data.medium;
    end
    function beam = set.medium(beam, val)
      beam.data.medium = val;
      beam.incident_beam.medium = val;
    end

    function power = get.power(beam)
      power = beam.total_beam.power;
    end
    function beam = set.power(beam, val)
      error('Changing scattered beam power not supported');
    end

    function tbeam = get.total_beam(beam)
      if strcmpi(beam.data_type, 'scattered')
        assert(~isempty(beam.incident_beam), ...
            'Need incident beam for conversion');
        tbeam = beam.data + beam.incident_beam;
      else strcmpi(beam.data_type, 'total')
        tbeam = beam.data;
      end
    end

    function tbeam = get.scattered_beam(beam)
      if strcmpi(beam.data_type, 'scattered')
        tbeam = beam.data;
      else strcmpi(beam.data_type, 'total')
        assert(~isempty(beam.incident_beam), ...
            'Need incident beam for conversion');
        tbeam = beam.data - beam.incident_beam;
      end
    end
  end
end
