classdef Scattered < ott.beam.abstract.Beam
% Represents a scattered beam.
% Inherits from :class:`abstract.Beam`.
%
% The scattered beam class has an additional property for the incident
% beam, which may be set to empty or a Abstract beam instance.
%
% Properties
%   - incident_beam   -- The incident ray object (can be set to [])
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

  properties
    incident_beam       % The incident ray object (can be set to [])
    type                % Type of beam ('scattered', 'total' or 'internal')
  end

  properties (Dependent)
    total_beam          % Instance of the beam with total type
    scattered_beam      % Instance of the beam with scattered type
  end

  methods
    function beam = Scattered(incident_beam, type, varargin)
      % Construct a new scattered beam representation
      %
      % Usage
      %   beam = Scattered(incident_beam, type, ...)
      %
      % Parameters
      %   - incident_beam ([]|Beam) -- Incident beam or emtpy.
      %
      %   - type (enum) -- Type of scattered beam.
      %     Either 'scattered', 'total' or 'internal'.
      %
      % Other parameters are passed to :class:`abstract.Beam`.

      beam = beam@ott.beam.abstract.Beam(varargin{:});
      beam.incident_beam = incident_beam;
      beam.type = type;
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
        tbeam.type = 'total';
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
        sbeam.type = 'scattered';
      else
        error('Unable to convert to specified type');
      end
    end
  end

  methods % Getters/setters
    function beam = set.incident_beam(beam, val)
      assert(isempty(val) || isa(val, 'ott.beam.abstract.Beam'), ...
        'Incident beam must be empty or a Ray object');
      beam.incident_beam = val;
    end

    function beam = set.type(beam, val)
      assert(any(strcmpi(val, {'scattered', 'total'})), ...
          'type must be ''scattered'', ''total'' or ''internal''');
      beam.type = val;
    end

    function tbeam = get.total_beam(beam)
      if strcmpi(beam.type, 'scattered')
        assert(~isempty(beam.incident_beam), ...
            'Need incident beam for conversion');
        tbeam = beam.scatteredField();
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
        tbeam = beam.totalField();
      else
        error('Unable to convert to specified type');
      end
    end
  end
end
