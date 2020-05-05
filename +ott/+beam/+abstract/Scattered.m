classdef Scattered < ott.beam.abstract.Abstract
% Represents a scattered beam.
% Inherits from :class:`Abstract`.
%
% The scattered beam class has an additional property for the incident
% beam, which may be set to empty or a Abstract beam instance.
%
% Properties
%   - incident_beam   -- The incident ray object (can be set to [])

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    incident_beam       % The incident ray object (can be set to [])
  end

  methods % Getters/setters
    function beam = set.incident_beam(beam, val)
      assert(isempty(val) || isa(val, 'ott.beam.abstract.Abstract'), ...
        'Incident beam must be empty or a Ray object');
      beam.incident_beam = val;
    end
  end
end
