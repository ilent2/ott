classdef Scattered < ott.beam.abstract.Beam
% Represents a scattered beam.
% Inherits from :class:`abstract.Beam`.
%
% The scattered beam class is a container which encapsulates a scattered
% beams type and the incident beam that cause the scattering.  Scattered
% beams behave like other beams, all properties are differed to the
% internal beam_data property.
%
% Properties
%   - beam_data       -- Beam object associated with the scattering
%   - incident_beam   -- The incident ray object (can be set to [])
%   - particle        -- Particle that caused the scattering (can be [])
%   - type            -- Type of beam ('scattered', 'total' or 'internal')
%
% Dependent properties
%   - total_beam      -- Instance of the beam with total type
%   - scattered_beam  -- Instance of the beam with scattered type
%
% Methods
%   - totalField      -- Calculate the total field beam
%   - scatteredField  -- Calculate the scattered field beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

% TODO: This class shouldn't inherit from abstract beam, since it has
%   no proeprties of its own.  It should implement the beam interface though.

  properties
    beam_data           % Beam object associated with the scattering
    incident_beam       % The incident ray object (can be set to [])
    particle            %  Particle that caused the scattering (can be [])
    type                % Type of beam ('scattered', 'total' or 'internal')
  end

  % TODO: Dependent beam properties or subsref/subsasgn methods

  properties (Dependent)
    total_beam          % Instance of the beam with total type
    scattered_beam      % Instance of the beam with scattered type
  end

  methods
    function beam = Scattered(type, beam, varargin)
      % Construct a new scattered beam representation
      %
      % Usage
      %   beam = Scattered(type, beam, ...)
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

      p = inputParser;
      p.addParameter('incident_beam', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.abstract.Beam(unmatched{:});
      beam.incident_beam = p.Results.incident_beam;
      beam = beam.setType(type);
    end

    function tbeam = totalField(beam, ibeam)
      % Calculate the total beam specifying the incident beam
      %
      % Usage
      %   total_beam = beam.totalField(incident_beam)
      %
      % If incident beam is not supplied, uses the internal incident_beam.
      %
      % If the beam is already a total beam, simply returns it
      % and raises a warning.

      if nargin == 1
        assert(~isempty(beam.incident_beam), ...
            'incident_beam must be specified in function or class');
        ibeam = beam.incident_beam;
      end

      assert(isa(ibeam, 'ott.beam.Beam'), ...
          'incident_beam must be an ott.beam.Beam');

      if strcmpi(beam.type, 'scattered')
        tbeam = 2*beam + ibeam;
        tbeam = tbeam.setType('total');
      elseif strcmpi(beam.type, 'total')
        tbeam = beam;
        if nargin == 2
          warning('Beam is already total');
        end
      else
        error('Unable to convert to specified type');
      end
    end

    function sbeam = scatteredField(beam, ibeam)
      % Calculate the scattered beam specifying the incident beam
      %
      % Usage
      %   scattered_beam = beam.scatteredField(incident_beam)
      %
      % If incident beam is not supplied, uses the internal incident_beam.
      %
      % If the beam is already a scattered beam, simply returns it
      % and raises a warning.

      if nargin == 1
        assert(~isempty(beam.incident_beam), ...
            'incident_beam must be specified in function or class');
        ibeam = beam.incident_beam;
      end

      assert(isa(ibeam, 'ott.beam.Beam'), ...
          'incident_beam must be an ott.beam.Beam');

      if strcmpi(beam.type, 'scattered')
        sbeam = beam;
        if nargin == 2
          warning('Beam is already scattered');
        end
      elseif strcmpi(beam.type, 'total')
        sbeam = 0.5*(beam - ibeam);
        sbeam = sbeam.setType('scattered');
      else
        error('Unable to convert to specified type');
      end
    end

    function beam = setType(beam, val)
      % Set the beam type paramter (without raising a warning)
      %
      % Usage
      %   beam = beam.setType(val);

      S = warning('off', 'ott:beam:abstract:Scattered:type_change');
      beam.type = val;
      warning(S);

      assert(nargout == 1, 'Too few outputs for setType');
    end
  end

  methods % Getters/setters
    function beam = set.incident_beam(beam, val)
      assert(isempty(val) || isa(val, 'ott.beam.abstract.Beam'), ...
        'Incident beam must be empty or a Ray object');
      beam.incident_beam = val;
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

    function tbeam = get.total_beam(beam)
      if strcmpi(beam.type, 'scattered')
        assert(~isempty(beam.incident_beam), ...
            'Need incident beam for conversion');
        tbeam = beam.totalField();
      elseif strcmpi(beam.type, 'total')
        tbeam = beam.beam_data;
      else
        error('Unable to convert to specified type');
      end
    end

    function tbeam = get.scattered_beam(beam)
      if strcmpi(beam.type, 'scattered')
        tbeam = beam.beam_data;
      elseif strcmpi(beam.type, 'total')
        assert(~isempty(beam.incident_beam), ...
            'Need incident beam for conversion');
        tbeam = beam.scatteredField();
      else
        error('Unable to convert to specified type');
      end
    end
  end
end
