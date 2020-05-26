classdef AbstractArray < ott.beam.properties.Beam
% Defines properties for abstract beam arrays.
% Inherits from :class:`Beam`.
%
% Properties
%   - beams         -- Internal array of abstract beams

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    beams           % Internal array of abstract beams
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
      beam.beams = val;
    end
  end
end
