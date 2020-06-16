classdef IntersectRayMarch
% Implement intersection methods using ray marching.
%
% Normals are calculated using the normalsXyz function.
%
% Abstract methods
%   - insideXyz
%   - normalsXyz
%   - intersectBoundingBox
%
% Methods
%   - intersectInternal
%   - intersectAllInternal

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Abstract)
    intersectBoundingBox
    normalsXyz
    insideXyz
  end

  methods (Hidden)
    function [locs, norms, dist] = intersectAllInternal(shape, vecs, varargin)
      % Calls intersectInternal repeatedly to find all intersections

      [locs, norms] = shape.intersectInternal(vecs);

      allLocs = reshape(locs, 3, 1, []);
      allNorms = reshape(norms, 3, 1, []);

      while ~all(any(isnan(allLocs(:, end, :)), 1), 3)

        % Get locs that aren't finished
        mask = ~any(isnan(squeeze(allLocs(:, end, :))), 1);
        ovecs = vecs.setData(vecs.origin(:, mask), vecs.direction(:, mask));

        [locs, norms] = shape.intersectInternal(ovecs);

        % Extend array
        allLocs(:, end+1, :) = nan;
        allNorms(:, end+1, :) = nan;

        % Assign new results
        allLocs(:, end, mask) = locs;
        allNorms(:, end, mask) = norms;

      end

      if nargout >= 3
        dist = locs - reshape(vecs.origin, 3, 1, []);
      end
    end

    function [locs, norms] = intersectInternal(shape, vecs, varargin)
      % Calculate intersection by ray marching

      % Minimum distance before intersection
      min_distance = p.Results.min_distance;

      dx = p.Results.stepSize;
      if isempty(dx)
        dx = shape.maxRadius ./ 100;
        if ~isfinite(dx)
          dx = 1.0e-3;
        end
      end

      % Start by calculating intersects with bounding box
      ints = shape.intersectBoundingBox(vecs);

      % Handle point inside (or on the edge)
      mask = isnan(ints(1, :));
      ints(1:3, mask) = vecs.origin(:, mask) ...
          + vecs.direction(:, mask) .* min_distance ...
          ./ vecnorm(vecs.direction(:, mask));

      % Calculate distance we want to search for each vector
      search_distance = vecnorm(ints(4:6, :) - ints(1:3, :));

      % Distance from bounding box to intersection
      locs = zeros(1, numel(vecs));
      orgs = ints(1:3, :);
      found = false(1, numel(vecs));
      dirs = vecs.direction ./ vecnorm(vecs.direction);

      % Determine which points were already inside
      was_inside = shape.insideXyz(orgs);

      % Find a point between the intersects in the shape
      remaining = locs < search_distance & ~found;
      while any(remaining)

        % March the ray
        locs(remaining) = locs(remaining) + dx;

        % Determine if new point is inside
        inside = shape.insideXyz(orgs(:, remaining) ...
          + locs(remaining).*dirs(:, remaining));
        found(remaining) = inside ~= was_inside(remaining);

        % Determine which points remain
        remaining = locs < search_distance & ~found;
      end

      % Remove logs that weren't found
      locs(~found) = nan;

      % Convert length to location vector
      locs = ints(1:3, :) + locs .* dirs;

      if nargout > 1
        norms = shape.normalsXyz(locs);
      end
    end
  end
end

