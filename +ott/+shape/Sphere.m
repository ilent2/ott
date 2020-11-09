classdef Sphere < ott.shape.Shape ...
    & ott.shape.mixin.AxisymStarShape ...
    & ott.shape.mixin.IsosurfSurfPoints ...
    & ott.shape.mixin.IntersectMinAll ...
    & ott.shape.mixin.IsSphereAbsProp
% Spherical particle.
% Inherits from :class:`Shape`.
%
% Properties
%   - radius        -- Radius of the sphere
%
% Additional properties inherited from base.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    xySymmetry         % True if the particle is xy-plane mirror symmetric
    perimeter          % Perimeter from in axis plane
    isSphere           % (Constant: true)
  end

  methods
    function shape = Sphere(varargin)
      % Construct a new sphere
      %
      % Usage
      %   shape = Sphere(radius, ...)
      %
      % Additional parameters passed to base class.

      p = inputParser;
      p.addOptional('radius', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shape.Shape(unmatched{:});
      radius = p.Results.radius;

      Nradius = numel(radius);
      assert(numel(shape) == 1 || Nradius == 1 ...
          || Nradius == numel(shape), ...
          'length of radius must match length of position');

      if numel(shape) == 1 && Nradius ~= 1
        shape = repelem(shape, 1, Nradius);
      end

      % Assign radius
      if Nradius > 1
        for ii = 1:Nradius
          shape(ii).radius = radius(ii);
        end
      else
        [shape.radius] = deal(radius);
      end
    end

    function shape = ott.shape.Ellipsoid(shape)
      % Convert object to a ellipsoid

      shape = ott.shape.Ellipsoid(shape.radius*[1,1,1]);
    end

    function shape = ott.shape.Superellipsoid(shape)
      % Convert object to a superellipsoid

      shape = ott.shape.Superellipsoid(shape.radius*[1,1,1], 1, 1);
    end

    function r = starRadii(shape, theta, phi)
      % Return the sphere radii
      %
      % Usage
      %   r = shape.starRadii(theta, phi)
      
      assert(all(size(theta) == size(phi)) ...
          || isscalar(theta) || isscalar(phi), ...
          'theta and phi must be same size or one must be scalar');
      sz = max(size(theta), size(phi));

      r = ones(sz) * shape.radius;
    end
  end

  methods (Hidden)
    function nxyz = normalsRtpInternal(~, rtp)
      % Calculate normals at the specified surface locations
      nxyz = ott.utils.rtpv2xyzv(repmat([1;0;0], 1, size(rtp, 2)), rtp);
    end

    function nxyz = normalsXyzInternal(~, xyz)
      nxyz = xyz ./ vecnorm(xyz);
    end

    function nxz = normalsTInternal(~, theta)
      nxz = [sin(theta); cos(theta)];
    end

    function [locs, norms, dist] = intersectAllInternal(shape, x0, x1)
      % Based on OTGO Spherical.m/intersectionpoint

      D = x1 - x0;
      Q1 = x0;

      A = dot(D, D);
      B = 2.*dot(Q1, D);
      C = dot(Q1, Q1) - shape.radius.^2;

      delta = B.^2 - 4*A.*C;

      t1 = (-B - sqrt(delta))./(2*A);
      t2 = (-B + sqrt(delta))./(2*A);

      % Remove rays that don't intersect
      t1(delta < 0) = nan;
      t2(delta < 0) = nan;

      % Remove rays in negative direction
      t1(t1 < 0) = nan;
      t2(t2 < 0) = nan;

      dist = cat(3, t1, t2);
      locs = Q1 + dist.*D;
      norms = shape.normalsXyz(reshape(locs, 3, []));

      % Ensure outputs have correct shape/size
      dist = permute(dist, [1, 3, 2]);
      locs = permute(locs, [1, 3, 2]);
      norms = permute(reshape(norms, 3, [], 2), [1, 3, 2]);

    end

    function shape = scaleInternal(shape, sc)
      shape.radius = shape.radius * sc;
    end
  end

  methods % Setters/getters
    function shape = set.radius(shape, val)
      assert(isnumeric(val) && isscalar(val), ...
          'radius must be numeric scalar');
      shape.radius = val;
    end

    function r = get.maxRadius(shape)
      % Calculate the maximum particle radius
      r = shape.radius;
    end

    function v = get.volume(shape)
      % Calculate the volume
      v = 4./3.*pi.*shape.radius.^3;
    end

    function bb = get.boundingBox(shape)
      bb = [-1, 1; -1, 1; -1, 1].*shape.radius;
    end

    function b = get.xySymmetry(~)
      b = true;
    end

    function p = get.perimeter(shape)
      % Calculate the perimiter of the object
      p = 2.0 * pi * shape.radius;
    end

    function b = get.isSphere(~)
      b = true;
    end
  end
end
