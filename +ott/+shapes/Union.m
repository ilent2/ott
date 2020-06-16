classdef Union < ott.shapes.Set
% Represents union between two shapes (| operator).
% Inherits from :class:`Set`.
%
% Properties
%   - shapes        -- Shapes forming the set
%   - maxRadius     -- Estimated from bounding box
%   - volume        -- Calculated numerically
%   - boundingBox   -- Surrounding all shapes
%
% Methods
%   - operator|     -- (Overloaded) Allow daisy chains

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    boundingBox         % Surrounds all shapes
  end

  methods
    function shape = Union(varargin)
      % Construct a union shape set.
      %
      % Usage
      %   shape = Union(shapes, ...)
      %
      % All parameters are passed to base class.

      shape = shape@ott.shapes.Set(varargin{:});
    end

    function sset = or(a, b)
      % Combine unions together smartly
      %
      % If position/rotation are different, creates a new union, otherwise
      % uses the existing union(s).
      % All other properties are copied from the first union.

      if isa(a, 'ott.shapes.Union') && isa(b, 'ott.shapes.Union') ...
          && a.hasMoved == false && b.hasMoved == false
        sset = a;
        sset.shapes = [sset.shapes, b.shapes];
      elseif isa(a, 'ott.shapes.Union') && a.hasMoved == false
        sset = a;
        sset.shapes = [sset.shapes, b];
      elseif isa(b, 'ott.shapes.Union') && b.hasMoved == false
        sset = b;
        sset.shapes = [a, sset.shapes];
      else
        sset = ott.shapes.Union([a, b]);
      end
    end
  end

  methods (Hidden)
    function b = insideRtpInternal(shape, rtp)
      % Determine if point is inside union of shapes

      b = shape.shapes(1).insideRtp(rtp, 'origin', 'global');
      for ii = 2:numel(shape.shapes)

        % We only need to check points that are not inside
        mask = ~b;

        b(mask) = b(mask) ...
            | shape.shapes(ii).insideRtp(rtp(:, mask), 'origin', 'global');
      end
    end

    function b = insideXyzInternal(shape, xyz)
      % Determine if point is inside union of shapes

      b = shape.shapes(1).insideXyz(xyz, 'origin', 'global');
      for ii = 2:numel(shape.shapes)
        % We only need to check points that are not inside
        mask = ~b;

        b(mask) = b(mask) ...
            | shape.shapes(ii).insideXyz(xyz(:, mask), 'origin', 'global');
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

      bb = [min(bbs, [], 2), max(bbs, [], 2)];
    end
  end
end
