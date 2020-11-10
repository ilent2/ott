classdef Intersection < ott.shape.Set
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
    maxRadius           % Radius surrounding all shapes
  end

  methods
    function shape = Intersection(varargin)
      % Construct a intersection shape set.
      %
      % Usage
      %   shape = Intersection(shapes, ...)
      %
      % All parameters are passed to base class.

      shape = shape@ott.shape.Set(varargin{:});
    end

    function sset = and(a, b)
      % Combine intersections together smartly
      %
      % If position/rotation are different, creates a new union, otherwise
      % uses the existing union(s).
      % All other properties are copied from the first union.

      if isa(a, 'ott.shape.Intersection') ...
          && isa(b, 'ott.shape.Intersection') ...
          && hasMoved(a) == false && hasMoved(b) == false
        sset = a;
        sset.shapes = [sset.shapes, b.shapes];
      elseif isa(a, 'ott.shape.Intersection') && hasMoved(a) == false
        sset = a;
        sset.shapes = [sset.shapes, b];
      elseif isa(b, 'ott.shape.Intersection') && hasMoved(b) == false
        sset = b;
        sset.shapes = [a, sset.shapes];
      else
        sset = ott.shape.Intersection([a, b]);
      end
      
      % Duplicated with ott.shape.Union, might move in future
      function b = hasMoved(shp)
        b = false;
        b = b & all(shp.position == [0;0;0]);
        b = b & all(all(shp.rotation == eye(3)));
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

    function r = get.maxRadius(shape)
      r = max(vecnorm(shape.boundingBox));
    end
  end
end
