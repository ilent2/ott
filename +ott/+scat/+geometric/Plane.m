classdef Plane < ott.scat.planewave.Plane
% Describes scattering of Rays by plane surfaces.
% Inherits from :class:ott.scat.planewave.Plane`.
%
% This class provides a similar scatter methods to its superclass
% except the returned type is rays and the position is updated
% to the intersection location.
%
% Methods
%   - scatter         -- Calculate scattered plane wave beams

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function [rbeam, tbeam] = scatter(plane, beam)
      % Calculate reflected and transmitted beams

      % Call the base method
      [rbeam, tbeam] = scatter@ott.scat.planewave.Plane(plane, beam);

      % Cast the result to rays
      rbeam = ott.beam.ScatteredRay('total', ...
          'like', rbeam, 'incident_beam', beam);
      tbeam = ott.beam.ScatteredRay('total', ...
          'like', tbeam, 'incident_beam', beam);

      % Calculate the intersection location
      locs = plane.intersect(beam);
      rbeam.origin = locs;
      tbeam.origin = locs;

    end
  end
end
