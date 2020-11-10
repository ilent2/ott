classdef Strata < ott.shape.Plane
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
      %   shape = Slab(depths, normal, ...)
      %
      % Optional named parameters
      %   - depths (numeric) -- Depth of each layer.
      %     Default: ``[0.2, 0.5]``.
      %
      %   - normal (3x1 numeric) -- Surface normal.
      %     Default: ``[]``.  Overwrite any value set with ``rotation``.
      %
      %   - offset (numeric) -- Offset of the plane from the position.
      %     Default: ``[]``.  Overwrites any values set with ``position``.
      %
      %   - position (3x1 numeric) -- Position of the plane.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Plane orientation.
      %     Default: ``eye(3)``.

      p = inputParser;
      p.addOptional('depths', [0.2, 0.5]);
      p.addOptional('normal', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shape.Plane('normal', p.Results.normal, unmatched{:});
      shape.depths = p.Results.depths;
    end

    function shape = ott.shape.Slab(oldshape)
      % Can be cast to a slab if the number of surfaces is 2

      assert(numel(oldshape.depths) == 1, 'Depth must have 1 element');

      shape = ott.shape.Slab('normal', oldshape.normal, ...
          'depth', oldshape.depths, ...
          'position', oldshape.offset .* oldshape.normal);
    end

    function shape = ott.shape.PatchMesh(shape, varargin)
      % Convert the plane into a patch
      %
      % Usage
      %   shape = ott.shape.PatchMesh(shape, ...)
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

      % Create list of verts
      verts = zeros(3, 4*(length(shape.depths)+1));
      verts(:, 1:4) = [v1 + v2, v1 - v2, -v1 - v2, -v1 + v2];
      for ii = 1:length(shape.depths)
        verts(:, (1:4) + ii*4) = verts(:, 1:4) ...
            + shape.depths(ii).*shape.normal;
      end
      
      % Create list of faces
      faces = reshape(1:4*(length(shape.depths)+1), 4, []);

      shape = ott.shape.PatchMesh(verts, faces, ...
          'position', shape.position, 'rotation', eye(3), unmatched{:});
    end
  end

  methods (Hidden)
    function [locs, norms, dist] = intersectAllInternal(shape, x0, x1)
      % Call plane method for each layer

      locs = zeros(3, length(shape.depths)+1, size(x0, 2));
      norms = zeros(3, length(shape.depths)+1, size(x0, 2));
      dist = zeros(1, length(shape.depths)+1, size(x0, 2));
      
      % First surface
      [locs(:, 1, :), norms(:, 1, :), dist(:, 1, :)] = ...
            intersectAllInternal@ott.shape.Plane(shape, x0, x1);

      for ii = 1:length(shape.depths)
        [locs(:, ii+1, :), norms(:, ii+1, :), dist(:, ii+1, :)] = ...
            intersectAllInternal@ott.shape.Plane(shape, ...
            x0 - shape.depths(ii)*shape.normal, ...
            x1 - shape.depths(ii)*shape.normal);
        locs(:, ii+1, :) = locs(:, ii+1, :) + shape.depths(ii)*shape.normal;
        dist(:, ii+1, :) = dist(:, ii+1, :) ...
          + vecnorm(shape.depths(ii)*shape.normal.*(x1 - x0))./vecnorm(x1 - x0);
      end
    end

    function shape = scaleInternal(shape, sc)
      shape.depths = shape.depths * sc;
      shape = scaleInternal@ott.shape.Plane(shape, sc);
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
