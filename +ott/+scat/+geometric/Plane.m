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

% using/distributing this file

  methods
    function plane = Plane(varargin)
      % Construct a new geometric plane scatterer
      %
      % Usage
      %   plane = Plane(normal, index_relative, ...)
      %
      % Arguments are passed to ott.scat.planewave.Plane constructor.

      plane = plane@ott.scat.planewave.Plane(varargin{:});
    end
  end

  methods (Hidden)
    function [rbeam, tbeam] = scatterInternal(plane, beam)
      % Calculate reflected and transmitted beams

      % Call the base method
      [rbeam, tbeam] = scatterInternal@ott.scat.planewave.Plane(plane, beam);

      % Cast the result to rays
      % TODO: This should just be a simple cast
      %   (once bug in planewave is fixed)
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
