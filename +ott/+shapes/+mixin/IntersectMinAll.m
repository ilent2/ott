classdef (InferiorClasses = {?ott.shapes.mixin.IntersectTriMesh}) ...
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
    function [locs, norms] = intersectInternal(shape, vecs, varargin)

      % Calculate both intersections and prune
      [locs, norms, dist] = shape.intersectAllInternal(vecs, varargin{:});

      % Find closest location
      [~, idx] = min(dist, [], 2);
      locs = locs(:, idx);
      norms = norms(:, idx);
    end
  end
end
