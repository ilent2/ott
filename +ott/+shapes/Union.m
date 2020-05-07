classdef Union < ott.shapes.Shape
% Represents union between two shapes.
% Inherits from :class:`+ott.+shapes.+Shape`.
%
% A point is considered to be inside the union if the point is inside
% any of the shapes in the union.
%
% Methods
%   inside    -- Determine if point is inside any contained shape.
%
% Properties
%   shapes    -- Shapes contained in this union
%   volume    -- Estimate of shape volume from sum of shapes in set
%   maxRadius -- Estimate of shape maximum radius
%
% See also Union

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% TODO: Normal calculation method

  properties
    shapes      % Shapes contained in this union
  end

  methods
    function shape = Union(shapes, varargin)
      % Construct a new union of shapes.
      %
      % Usage
      %   shape = Union(shape_array, ...)
      %
      % Parameters
      %   shape_array -- Array of :class:`ott.shapes.Shapes`.
      %
      % Optional named arguments
      %   - position (3 numeric) -- Position of the shape.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Orientation of the shape.
      %     Default: ``eye(3)``.

      shape = shape@ott.shapes.Shape(varargin{:});

      % Store shapes
      [shape.shapes] = deal(shapes);
    end
  end

  methods (Hidden)

    function b = insideRtpInternal(shape, rtp, varargin)
      % Determine if point is inside union of shapes

      b = shape.shapes(1).insideRtpInternal(rtp, 'origin', 'world');
      for ii = 2:numel(shape.shapes)
        % TODO: We only need to check points that are not inside
        b = b | shape.shapes(ii).insideRtpInternal(rtp, 'origin', 'world');
      end
    end

    function b = insideXyzInternal(shape, xyz, varargin)
      % Determine if point is inside union of shapes

      b = shape.shapes(1).insideXyzInternal(rtp, 'origin', 'world');
      for ii = 2:numel(shape.shapes)
        % TODO: We only need to check points that are not inside
        b = b | shape.shapes(ii).insideXyzInternal(rtp, 'origin', 'world');
      end
    end

    function radius = get_maxRadius(shape, varargin)
      % Calculate maximum distance from shape centre
      %
      % Estimates the maximum radius based on the distance of each
      % shape from the origin and each shapes radius.

      radius = vecnorm(shape.shapes(1).position, 2) ...
          + shape.shapes(1).maxRadius;
      for ii = 2:numel(shape.shapes)
        radius = max(radius, vecnorm(shape.shapes(ii).position, 2) ...
          + shape.shapes(ii).maxRadius);
      end
    end

    function volume = get_volume(shape, varargin)
      % Calculate the volume of the combined shape
      %
      % Estimates the volume as the sum of the volume of each
      % child shape.

      volume = shape.shapes(1).volume;
      for ii = 2:numel(shape.shapes)
        volume = volume + shape.shapes(ii).volume;
      end
    end
  end

  methods % Getters/setters
    function shape = set.shapes(shape, value)
      assert(numel(value) >= 1, 'Number of shapes must be >= 1');
      assert(all(isa(value, 'ott.shapes.Shape')), ...
          'Shapes must all inherit from ott.shapes.Shape');
      shape.shapes = value;
    end
  end
end
