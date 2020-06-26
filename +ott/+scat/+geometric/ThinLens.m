classdef ThinLens
% Thin lens approximation.
% Inherits from :class:`ott.shapes.Plane`.
%
% The thin lens focuses rays to a focal plane parallel to the lens plane.
%
% For a ray incident normal on the lens, the ray is focussed to a
% point a distance `focal_length` away from the lens centre.
% For non-normal rays, the focal spot shifts depending on the
% angle of the ray relative to the lens.
%
% Properties
%   - normal            -- Normal to the lens surface
%   - position          -- Position of the lens
%   - focal_length      -- Focal length of the lens (can be negative)
%
% Methods
%   - intersection      -- Calculates the intersection locations
%   - scatter           -- Calculate the transmitted beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

% TODO: getGeometry method for lens surface

  properties
    focal_length       % Focal length of the lens (can be negative)
  end

  methods
    function lens = ThinLens(focal_length, varargin)
      % Construct a new thin lens approximation
      %
      % Usage
      %   lens = ThinLens(focal_length, ...)
      %
      % Parameters
      %   - focal_length (numeric) -- Focal length of the lens, can be
      %     negative.
      %
      % Optional named arguments
      %   - normal (3x1 numeric) -- Normal to the lens surface.
      %     Default: [0;0;1]
      %
      %   - position (3x1 numeric) -- Position of the lens.
      %     Default: [0;0;0]

      p = inputParser;
      p.addParameter('position', [0;0;0]);
      p.addParameter('normal', [0;0;1]);
      p.parse(varargin{:});

      lens = lens@ott.shapes.Plane(p.Results.normal, 0.0);
      lens.position = p.Results.position;
      lens.focal_length = focal_length;
    end

    function tbeam = scatter(lens, ibeam)
      % Calculate the transmitted beam
      %
      % Usage
      %   tbeam = lens.scatter(ibeam)
      %
      % Calculates the beam transmitted through the thin lens.

      % Calculate intersection locations
      locs = lens.intersect(ibeam);

      % Calculate focal point location
      dirs = ibeam.direction ./ vecnorm(ibeam.direction);
      ndirs = sum(dirs .* lens.normal, 1);
      fp = lens.position + lens.focal_length .* dirs ./ ndirs;

      % Calculate new directions
      directions = sign(lens.focal_length)*(fp - locs);
      directions = directions ./ vecnorm(directions);

      % Create transmitted beam
      tbeam = ott.beam.ScatteredRay(ibeam, 'origin', locs, ...
          'direction', directions, 'like', ibeam);
    end
  end

  methods % Getters/setters
    function lens = set.focal_length(lens, val)
      assert(isnumeric(val) && isscalar(val), ...
        'focal length must be numeric scalar');
      lens.focal_length = val;
    end
  end
end

