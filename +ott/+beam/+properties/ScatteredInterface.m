classdef ScatteredInterface < ott.beam.properties.Beam
% Interface for scattered beam properties
% Inherits from :class:`Beam`.
%
% Abstract properties
%   - incident_beam -- Incident beam
%   - total_beam      -- Total-field beam (incident+scattered)
%   - scattered_beam  -- Scattered-field beam
%   - type          -- Default type used for casts/visualisations
%   - omega           -- omega property of beam_data
%   - medium          -- medium property of beam_data
%   - power           -- Total beam power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    type
  end

  properties (Abstract)
    incident_beam
    total_beam
    scattered_beam
    omega               % omega property of beam_data
    medium              % medium property of beam_data
    power               % Total beam power
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Adds similar properties to an array
      if isa(other, 'ott.beam.properties.ScatteredInterface')
        args = ott.utils.addDefaultParameter(...
            'type', other.type, args);
      end
      args = ott.beam.properties.Beam.likeProperties(...
          other, args);
    end
  end

  methods
    function beam = ScatteredInterface(varargin)
      % Construct scattered interface properties
      %
      % Usage
      %   ScatteredInterface(type, ...)

      p = inputParser;
      p.addOptional('type', [], ...
          @(x) any(strcmpi(x, {'scattered', 'total'})));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Beam(unmatched{:});
      beam.type = p.Results.type;
    end
  end

  methods % Getters/setters
    function beam = set.type(beam, val)
      assert(any(strcmpi(val, {'scattered', 'total'})), ...
          'type must be ''scattered'' or ''total''');
      beam.type = val;
    end
  end
end
