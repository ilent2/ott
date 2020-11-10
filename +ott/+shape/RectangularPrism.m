classdef RectangularPrism < ott.shape.Shape ...
    & ott.shape.mixin.CoordsCart ...
    & ott.shape.mixin.Patch ...
    & ott.shape.mixin.IntersectMinAll
% Simple geometric rectangular prism.
% Inherits from :class:`Shape`.
%
% Properties
%   - widths        -- Widths of each side [x; y; z]
%
% Supported casts
%   - ott.shape.Cube -- Only if widths all match
%
% Additional properties inherited from base.

% Copyright 2018-2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    widths         % Widths of each side
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    starShaped         % True if the particle is star-shaped
    xySymmetry         % True if the particle is xy-plane mirror symmetric
    zRotSymmetry       % z-axis rotational symmetry of particle
  end

  properties (Dependent, SetAccess=protected)
    verts              % Vertices forming surface
    faces              % Faces forming surface
  end

  methods
    function shape = RectangularPrism(varargin)
      % Construct a rectangular base prism
      %
      % Usage
      %   shape = RectangularPrism(widths, ...)
      %   Parameters can be passed as named arguments.
      %
      % Additional parameters are passed to base.

      p = inputParser;
      p.addOptional('widths', [1.0, 2.0, 3.0]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shape.Shape(unmatched{:});
      shape.widths = p.Results.widths;
    end

    function shape = ott.shape.Cube(shape)
      % Cast shape to cube (only supported if all widths match

      assert(all(shape.widths(1) == shape.widths), ...
          'All widths must match for cast to cube');

      shape = ott.shape.Cube(shape.widths(1), ...
          'position', shape.position, ...
          'rotation', shape.rotation);
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz)
      assert(size(xyz, 1) == 3, 'xyz must be 3xN matrix');
      b = all(abs(xyz) <= shape.widths./2, 1);
    end

    function n = normalsXyzInternal(shape, xyz)
      % Defers to the cube normals method after scaling coordinates

      assert(size(xyz, 1) == 3, 'xyz must be 3xN matrix');
      xyz = xyz ./ shape.widths(:);
      n = ott.shape.Cube(1.0).normalsXyzInternal(xyz);
    end

    function [locs, norms, dist] = intersectAllInternal(shape, x0, x1)
      % Defer to cube after scaling coordinates

      cube = ott.shape.Cube(1.0);
      [locs, norms] = cube.intersectAllInternal(...
          x0 ./ shape.widths, x1 ./ shape.widths);

      % Rescale locations (normals are fine)
      locs = locs .* shape.widths;

      % Compute distances
      O = reshape(x0, 3, 1, []);
      D = reshape(x1 - x0, 3, 1, []);
      dist = sum((locs - O).*D, 1)./vecnorm(D);
    end

    function varargout = intersectInternal(shape, varargin)
      % Disambiguate
      [varargout{1:nargout}] = intersectInternal@ ...
          ott.shape.mixin.IntersectMinAll(shape, varargin{:});
    end

    function shape = scaleInternal(shape, sc)
      shape.widths = shape.widths * sc;
    end
  end

  methods % Getters/setters
    function shape = set.widths(shape, val)
      assert(isnumeric(val) && numel(val) == 3 && all(val >= 0), ...
          'widths must be positive numeric 3-vector');
      shape.widths = val(:);
    end

    function r = get.maxRadius(shape)
      r = vecnorm(shape.widths./2);
    end
    function shape = set.maxRadius(shape, ~) %#ok<INUSD>
      error('ott:shape:RectangularPrism:set_maxradius', ...
          'Cannot set maxRadius, set widths instead');
    end

    function v = get.volume(shape)
      v = prod(shape.widths);
    end
    function shape = set.volume(shape, ~) %#ok<INUSD>
      error('ott:shape:RectangularPrism:set_volume', ...
          'Cannot set volume, set widths instead');
    end

    function bb = get.boundingBox(shape)
      bb = [-1, 1; -1, 1; -1, 1].*shape.widths./2;
    end

    function b = get.starShaped(~)
      b = true;
    end
    function b = get.xySymmetry(~)
      b = true;
    end
    function q = get.zRotSymmetry(shape)
      if shape.widths(1) == shape.widths(2)
        q = 4;
      else
        q = 2;
      end
    end

    function verts = get.verts(shape)
      X = [1, 1, -1, -1, 1, 1, -1, -1];
      Y = [-1, 1, 1, -1, -1, 1, 1, -1];
      Z = [-1, -1, -1, -1, 1, 1, 1, 1];

      verts = [X(:), Y(:), Z(:)].' .* shape.widths/2;
    end

    function faces = get.faces(~)
      faces = [1, 2, 6, 5; 2, 3, 7, 6; 3, 4, 8, 7; 4, 1, 5, 8; ...
          2, 1, 4, 3; 5, 6, 7, 8].';
    end
  end
end
