classdef Strata < ott.shapes.Plane
% Shape describing a series of stratified interfaces.
% Inherits from :class:`Plane`.
%
% This shape describes a series of layered planes.  When the number
% of layers is equal to 2, this object can be converted to a Slab.
% All points above the first layer are considered to be inside the shape.
%
% Properties
%   - normal      -- Vector representing surface normal
%   - offset      -- Offset of surface from coordinate origin
%   - depths      -- Depth of each layer (must be positive)

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    depths        % Depth of each layer (must be positive)
  end

  methods
    function shape = Strata(varargin)
      % Construct a new infinite slab
      %
      % Usage
      %   shape = Slab(normal, depth)
      %
      %   shape = Slab(normal, offset, depth)
      %
      % Parameters
      %   - normal (3x1 numeric) -- Surface normal
      %   - depth (numeric) -- Depth of each layer.
      %   - offset (numeric) -- Offset of first surface to origin.

      p = inputParser;
      p.addOptional('normal', []);
      p.addOptional('depths', [0.2, 0.5]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shapes.Plane('normal', p.Results.normal, unmatched{:});
      shape.depths = p.Results.depths;
    end

    function shape = ott.shapes.Slab(oldshape)
      % Can be cast to a slab if the number of surfaces is 2

      assert(numel(oldshape.depths) == 1, 'Depth must have 1 element');

      shape = ott.shapes.Slab(oldshape.normal, oldshape.depth, ...
          'offset', p.Results.offset);
    end

    function shape = ott.shapes.PatchMesh(shape, varargin)
      % Convert the plane into a patch
      %
      % Usage
      %   shape = ott.shapes.PatchMesh(shape, ...)
      %
      % Optional named parameters
      %   - scale (numeric) -- Size of the generated patch.
      %
      % Unmatched parameters are passed to patch constructor.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      v1 = shape.rotation(:, 1);
      v2 = shape.rotation(:, 2);

      v1 = v1 ./ vecnorm(v1) .* p.Results.scale ./ 2;
      v2 = v2 ./ vecnorm(v2) .* p.Results.scale ./ 2;

      verts = [v1 + v2, v1 - v2, -v1 - v2, -v1 + v2];
      faces = [1; 2; 3; 4];

      % Duplicate for additional layers
      for ii = 1:length(shape.depths)
        verts = [verts, verts(:, 1:4) + shape.depths(ii).*shape.normal];
        faces = [faces, (1:4).'+ii*4];
      end

      shape = ott.shapes.PatchMesh(verts, faces, ...
          'position', shape.position, 'rotation', eye(3), unmatched{:});
    end
  end

  methods (Hidden)
    function [locs, norms, dist] = intersectAllInternal(shape, vecs)
      % Call plane method for each layer

      locs = zeros(3, length(shape.depths), length(vecs));
      norms = zeros(3, length(shape.depths), length(vecs));
      dist = zeros(1, length(shape.depths), length(vecs));

      for ii = 1:length(shape.depths)
        tvecs = ott.utils.Vector('origin', vecs.origin - shape.depths(ii), ...
            'direction', vecs.direction);
        [locs(:, ii, :), norms(:, ii, :), dist(:, ii, :)] = ...
            intersectAllInternal@ott.shapes.Plane(shape, tvecs);
      end
    end

    function shape = scaleInternal(shape, sc)
      shape.depths = shape.depths * sc;
      shape = scaleInternal@ott.shapes.Plane(sc);
    end
  end

  methods % Getters/setters
    function shape = set.depths(shape, val)
      assert(isvector(val) && all(val > 0), ...
        'depths must be vector of positive values');
      shape.depths = val;
    end
  end
end
