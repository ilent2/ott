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
    function beam = Scattered(type, varargin)
      % Construct a new scattered beam representation
      %
      % Usage
      %   beam = Scattered(type, ...)
      %
      % Parameters
      %   - incident_beam ([]|Beam) -- Incident beam or emtpy.
      %     Default: ``[]``.
      %
      %   - type (enum) -- Type of scattered beam.
      %     Either 'scattered', 'total' or 'internal'.
      %
      % Other parameters are passed to :class:`abstract.Beam`.
      
      p = inputParser;
      p.addParameter('incident_beam', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.abstract.Beam(unmatched{:});
      beam.incident_beam = p.Results.incident_beam;
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
    
    function beam = setType(beam, val)
      % Set the beam type paramter (without raising a warning)
      %
      % Usage
      %   beam = beam.setType(val);
      
      S = warning('off', 'ott:beam:abstract:Scattered:type_change');
      beam.type = val;
      warning(S);
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
        tbeam = beam.scatteredField();
      else
        error('Unable to convert to specified type');
      end
    end
  end
end
