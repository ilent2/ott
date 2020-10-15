classdef (Abstract) Set < ott.shape.Shape ...
    & ott.shape.mixin.NumericalVolume ...
    & ott.shape.mixin.VarXySymmetry ...
    & ott.shape.mixin.VarZRotSymmetry ...
    & ott.shape.mixin.VarStarShaped ...
    & ott.shape.mixin.IsosurfSurf ...
    & ott.shape.mixin.IntersectRayMarch
% Collection of shapes
% Inherits from :class:`Shape`.
%
% This is the base class for collections of shapes including
% :class:`Union` and :class:`Intersection`.
%
% Properties
%   - shapes        -- Shapes forming the set
%   - volume        -- Calculated numerically
%   - starShaped    -- Variable, default false
%   - xySymmetry    -- Variable, default false
%   - zRotSymmetry  -- Variable, default 1

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    shapes              % Shapes contained in the set
  end

  properties (Dependent)
    maxRadius           % Radius surrounding all shapes
  end

  methods
    function shape = Set(shapes, varargin)
      % Construct a new shape set
      %
      % This is the abstract constructor for shape sets.  Use Union
      % or Intersection for instances.
      %
      % Usage
      %   shape = shape@ott.shape.Set(shapes, ...)
      %
      % Stores shapes and passes optional arguments to base.

      shape = shape@ott.shape.Shape(varargin{:});
      shape.shapes = shapes;
    end
  end

  methods (Hidden)
    function shape = scaleInternal(shape, sc)
      shape.shapes = shape.shapes.scale(sc);
    end
  end

  methods % Getters/setters
    function shape = set.shapes(shape, val)
      assert(numel(val) >= 1, 'number of shapes must be >= 1');
      assert(isa(val, 'ott.shape.Shape'), ...
          'shapes must inherit from ott.shape.Shape');
      shape.shapes = val;
    end

    function r = get.maxRadius(shape)
      r = max(vecnorm(shape.boundingBox));
    end
  end
end

