classdef (Abstract) ScatteredStem < ott.beam.abstract.CastBoth
% Abstract base class (stem) for scattered beams
% Inherits from :class:`CastBoth` and
% :class:`ott.beam.properties.ScatteredInterface`.
%
% Casts
%   Beam      -- Retrieves the total_beam or scattered_beam

  properties (Abstract)
    type
    scattered_beam
    total_beam
  end

  methods
    function beam = ott.beam.Beam(beam, varargin)
      % Get total or scattered beam

      assert(isa(beam, 'ott.beam.abstract.ScatteredStem'), ...
        'First argument must be a beam');

      if strcmpi(beam.type, 'scattered')
        beam = beam.scattered_beam;
      else
        beam = beam.total_beam;
      end
    end
  end
end
