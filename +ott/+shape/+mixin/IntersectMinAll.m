classdef (InferiorClasses = {?ott.shape.mixin.IntersectTriMesh}) ...
    IntersectMinAll
% Calculate intersection by first calculating all intersections.
%
% Inferior classes: :class:`IntersectTriMesh`.
%
% Methods
%   - intersectInternal
%
% Abstract methods
%   - intersectAllInternal

  methods (Abstract)
    intersectAllInternal(obj)
  end

  methods (Hidden)
    function [locs, norms] = intersectInternal(shape, varargin)

      % Calculate both intersections and prune
      [locs, norms, dist] = shape.intersectAllInternal(varargin{:});

      % Find closest location
      [~, idx] = min(dist, [], 2);
      idx = sub2ind([size(locs, 2), size(locs, 3)], idx(:).', 1:numel(idx));
      locs = locs(:, idx);
      norms = norms(:, idx);
    end
  end
end
