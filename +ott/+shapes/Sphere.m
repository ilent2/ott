classdef Sphere < ott.shapes.Shape ...
    & ott.shapes.mixin.AxisymStarShape ...
    & ott.shapes.mixin.IsosurfSurfPoints ...
    & ott.shapes.mixin.IntersectMinAll
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

      shape = shape@ott.shapes.Shape(unmatched{:});
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

    function r = starRadii(shape, theta, phi)
      % Return the sphere radii
      %
      % Usage
      %   r = shape.starRadii(theta, phi)

      r = ones(size(theta)) * shape.radius;
    end

    function [rtp, n, ds] = boundarypoints(shape, varargin)
      % BOUNDARYPOINTS calculates boundary points for surface integral
      %
      % [rtp, n, ds] = BOUDNARYPOINTS(npts) calculates the boundary points
      % and surface normal vectors in spherical coordinates and the area
      % elements of each ring.
      %
      % BOUNDARYPOINTS('Nmax', Nmax) takes a guess at a suitable npts
      % for the given Nmax.

      ntheta = shape.boundarypoints_npts(varargin{:});

      % Sphere has equally spaced angles
      theta = 0.0:(pi/(ntheta-1)):pi;
      phi = zeros(size(theta));

      xyz = shape.locations(theta, phi);
      nxyz = xyz ./ shape.radius;

      % Convert cylindrical coordinates into spherical coordinates
      [n, rtp] = ott.utils.xyzv2rtpv(nxyz, xyz);

      % Calculate area elements
      ds = shape.boundarypoints_area(xyz(:, 1), xyz(:, 3), ...
          xyz(:, 1), xyz(:, 3), rtp);

    end
  end

  methods (Hidden)
    function nxyz = normalsRtpInternal(shape, rtp)
      % Calculate normals at the specified surface locations
      nxyz = ott.utils.rtpv2xyzv(repmat([1;0;0], 1, size(rtp, 2)), rtp);
    end

    function nxyz = normalsXyzInternal(shape, xyz)
      nxyz = xyz ./ vecnorm(xyz);
    end

    function nxz = normalsTInternal(shape, theta)
      nxz = [sin(theta); cos(theta)];
    end

    function [locs, norms, dist] = intersectAllInternal(shape, vecs)
      % Based on OTGO Spherical.m/intersectionpoint

      D = vecs.direction;
      Q1 = vecs.origin;

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

    function b = get.xySymmetry(shape)
      b = true;
    end

    function p = get.perimeter(shape)
      % Calculate the perimiter of the object
      p = 2.0 * pi * shape.radius;
    end
  end
end
