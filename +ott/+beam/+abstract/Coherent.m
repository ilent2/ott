classdef Coherent < ott.beam.properties.Coherent ...
    & ott.beam.abstract.Array
% Specialisation for coherent abstract beam arrays.
% Inherits from :class:`ott.beam.properties.Coherent` and :class:`Array`.
%
% Supported casts
%   - Coherent

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = ott.beam.Coherent(beam, varargin)
      % Cast the abstract array to a full array
      %
      % If the beam is an array of coherent arrays, casts to a
      % Coherent array of coherent beams.

      if numel(beam) > 1
        oldbeam = beam;
        beam = ott.beam.Coherent(size(beam));
        for ii = 1:numel(oldbeam)
          beam(ii) = ott.beam.Coherent(oldbeam(ii));
        end
      else
        beam = ott.beam.Coherent(beam.beams);
      end
    end
  end
end
