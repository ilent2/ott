classdef AbstractArray < ott.beam.properties.Beam
% Defines properties for abstract beam arrays.
% Inherits from :class:`Beam`.
%
% Properties
%   - beams         -- Internal array of abstract beams
%   - omega         -- Beam optical frequency
%   - medium        -- Beam optical medium
%   - power         -- Total beam power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    beams           % Internal array of abstract beams
  end

  properties (Dependent)
    omega     % Beam optical frequency
    medium    % Beam optical medium
  end

  methods
    function beam = AbstractArray(sz, varargin)
      % Construct a new beam array
      %
      % Usage
      %   beam = beam@ott.beam.properties.AbstractArray(sz, ...)
      %   Additional properties are passed to properties.Beam.

      beam = beam@ott.beam.properties.Beam(varargin{:});
      beam.beams = repelem(ott.beam.abstract.Empty, sz);
    end
  end

  methods % Getters/setters
    function beam = set.beams(beam, val)
      assert(isa(val, 'ott.beam.abstract.Beam'), ...
          'beams must be a abstract.Beam array');

      assert(numel(unique([val.omega])) == 1, ...
          'beams must all have same optical frequency');
      assert(numel(unique([val.medium])) == 1, ...
          'beams must all have same medium');

      beam.beams = val;
    end

    function omega = get.omega(beam)
      omega = beam.beams(1).omega;
    end
    function beam = set.omega(beam, val)
      [beam.beams.medium] = deal(val);
    end

    function medium = get.medium(beam)
      medium = beam.beams(1).medium;
    end
    function beam = set.medium(beam, val)
      [beam.beams.medium] = deal(val);
    end
  end
end
