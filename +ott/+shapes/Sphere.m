classdef Sphere < ott.shapes.Shape ...
    & ott.shapes.mixin.AxisymStarShape ...
    & ott.shapes.mixin.IsosurfSurfPoints
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
    function shape = Sphere(radius, varargin)
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
