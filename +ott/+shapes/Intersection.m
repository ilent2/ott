classdef Intersection < ott.shapes.Set
% Represents intersection between two shapes (& operator).
% Inherits from :class:`Set`.
%
% Properties
%   - shapes        -- Shapes forming the set
%   - maxRadius     -- Estimated from bounding box
%   - volume        -- Calculated numerically
%   - boundingBox   -- Surrounding intersection
%
% Methods
%   - operator&     -- (Overloaded) Allow daisy chains

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    boundingBox         % Surrounding intersection
  end

  methods
    function shape = Intersection(varargin)
      % Construct a intersection shape set.
      %
      % Usage
      %   shape = Intersection(shapes, ...)
      %
      % All parameters are passed to base class.

      shape = shape@ott.shapes.Set(varargin{:});
    end

    function sset = and(a, b)
      % Combine intersections together smartly
      %
      % If position/rotation are different, creates a new union, otherwise
      % uses the existing union(s).
      % All other properties are copied from the first union.

      if isa(a, 'ott.shapes.Intersection') ...
          && isa(b, 'ott.shapes.Intersection') ...
          && a.hasMoved == false && b.hasMoved == false
        sset = a;
        sset.shapes = [sset.shapes, b.shapes];
      elseif isa(a, 'ott.shapes.Intersection') && a.hasMoved == false
        sset = a;
        sset.shapes = [sset.shapes, b];
      elseif isa(b, 'ott.shapes.Intersection') && b.hasMoved == false
        sset = b;
        sset.shapes = [a, sset.shapes];
      else
        sset = ott.shapes.Intersection([a, b]);
      end
    end
  end

  methods (Hidden)
    function b = insideRtpInternal(shape, rtp)
      % Determine if point is inside union of shapes

      b = shape.shapes(1).insideRtp(rtp, 'origin', 'global');
      for ii = 2:numel(shape.shapes)
        b = b & shape.shapes(ii).insideRtp(rtp, 'origin', 'global');
      end
    end

    function b = insideXyzInternal(shape, xyz)
      % Determine if point is inside union of shapes

      b = shape.shapes(1).insideXyz(xyz, 'origin', 'global');
      for ii = 2:numel(shape.shapes)
        b = b & shape.shapes(ii).insideXyz(xyz, 'origin', 'global');
      end
    end
  end

  methods % Getters/setters
    function bb = get.boundingBox(shape)

      bbs = zeros(3, 2*numel(shape.shapes));
      for ii = 1:numel(shape.shapes)
        obb = shape.shapes(ii).getBoundingBox('origin', 'global');
        bbs(:, (1:2)+(ii-1)*2) = obb;
      end

      bb = [max(bbs(:, 1:2:end), [], 2), min(bbs(:, 2:2:end), [], 2)];
    end
  end
end
