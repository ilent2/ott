classdef Cylinder < ott.shape.Shape ...
    & ott.shape.mixin.AxisymStarShape ...
    & ott.shape.mixin.IsosurfSurfPoints ...
    & ott.shape.mixin.IntersectMinAll
% A simple cylindrical shape.
% Inherits from :class:`Shape`, :class:`mixin.AxisymStarShape` and
% :class:`mixin.IsosurfSurfPoints`.
%
% Properties
%  - radius       -- Radius of the cylinder
%  - height       -- Height of the cylinder

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius        % Radius of the cylinder
    height        % Height of the cylinder
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    xySymmetry         % True if the particle is xy-plane mirror symmetric
    perimeter          % Perimeter from in axis plane
  end

  methods
    function shape = Cylinder(varargin)
      % Construct a new cylinder
      %
      % Usage
      %   shape = Cylinder(radius, height, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Additional parameters passed to base :class:`Shape`.

      p = inputParser;
      p.addOptional('radius', 1.0);
      p.addOptional('height', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shape.Shape(unmatched{:});
      shape.radius = p.Results.radius;
      shape.height = p.Results.height;
    end

    function shape = ott.shape.PatchMesh(shape, varargin)
      % Cast shape to a PatchMesh
      %
      % Usage
      %   shape = ott.shape.PatchMesh(shape, ...)
      %
      % Optional named parameters
      %   - resolution ([ntheta, nphi]) -- Number of faces in theta/phi
      %     directions.  Default: ``[0, 20]`` uses cylinder corners for
      %     theta direction unless otherwise set.
      %
      % Additional named parameters are passed to constructor.

      p = inputParser;
      p.addParameter('resolution', [0, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      res = p.Results.resolution;
      if res(1) == 0
        edge_angle = atan(shape.height/(2*shape.radius));
        theta = [0, pi/2 - edge_angle, pi/2 + edge_angle, pi];
      else
        theta = linspace(0, 2*pi, res(1)+1);
      end
      phi = linspace(0, 2*pi, res(2)+1);
      [T, P] = meshgrid(theta, phi);

      % Calculate radii
      R = shape.starRadii(T(1:end-1, :), P(1:end-1, :));

      % Ensure ends and edge locations match
      % Reduces work and makes FromSurfMatrix behave nicer
      R(:, 1) = R(1, 1);
      R(:, end) = R(1, end);
      T(end, :) = T(1, :);
      P(end, :) = P(1, :);
      R(end+1, :) = R(1, :);

      % Convert to Cartesian coordinates
      [X, Y, Z] = ott.utils.rtp2xyz(R, T, P);

      shape = ott.shape.PatchMesh.FromSurfMatrix(X, Y, Z, ...
          'position', shape.position, 'rotation', shape.rotation);
    end

    function r = starRadii(shape, theta, phi)
      % Return the radii of the requested points
      %
      % Usage
      %   r = shape.starRadii(theta, phi)

      % Find angle of edges from xy plane
      edge_angle = atan(shape.height/(2*shape.radius));

      sides = ( theta >= pi/2 - edge_angle ) ...
          & ( theta <= pi/2 + edge_angle ) ;
      top = theta < pi/2 - edge_angle;
      bottom = theta > pi/2 + edge_angle;
      ends = top | bottom;

      r = zeros(size(theta));
      r(sides) = shape.radius ./ sin(theta(sides));
      r(ends) = (shape.height/2) ./ abs(cos(theta(ends)));
    end
  end

  methods (Hidden)
    function nxz = normalsTInternal(shape, theta)

      % Find angle of edges from xy plane
      edge_angle = atan(shape.height/(2*shape.radius));

      sides = ( theta >= pi/2 - edge_angle ) ...
          & ( theta <= pi/2 + edge_angle ) ;
      top = theta < pi/2 - edge_angle;
      bottom = theta > pi/2 + edge_angle;
      ends = top | bottom;

      nxz = zeros(2, numel(theta));
      nxz(1, sides) = sin(theta(sides));
      nxz(2, sides) = cos(theta(sides));
      nxz(1, ends) = abs(cos(theta(ends)));
      nxz(2, ends) = - sin(theta(ends)) .* sign(cos(theta(ends)));
    end

    function [locs, norms, dist] = intersectAllInternal(shape, x0, x1, varargin)
      % Based on OTGO Cylindrical.m/intersectionpoint

      Q1 = x0;
      D = x1 - x0;

      % Calculate intersection locations with cylindrical part

      A = dot(D(1:2, :), D(1:2, :));
      B = 2.*dot(Q1(1:2, :), D(1:2, :));
      C = dot(Q1(1:2, :), Q1(1:2, :)) - shape.radius.^2;
      delta = B.^2 - 4*A.*C;

      t1 = (-B - sqrt(delta))./(2*A);
      t2 = (-B + sqrt(delta))./(2*A);

      % Remove rays before origin
      t1(t1 < 0) = nan;
      t2(t2 < 0) = nan;

      P1 = x0 + t1.*(x1 - x0);
      P2 = x0 + t2.*(x1 - x0);

      % Remove rays that don't intersect
      P1(:, delta < 0) = nan;
      P2(:, delta < 0) = nan;

      % Calculate intersection locations with end caps
      t3 = (shape.height./2 - Q1(3, :))./D(3, :);
      t4 = (-shape.height./2 - Q1(3, :))./D(3, :);

      % Remove rays before origin
      t3(t3 < 0) = nan;
      t4(t4 < 0) = nan;

      P3 = x0 + t3.*(x1 - x0);
      P4 = x0 + t4.*(x1 - x0);

      % Remove rays that don't intersect
      P3(:, vecnorm(P3(1:2, :)) > shape.radius) = nan;
      P4(:, vecnorm(P4(1:2, :)) > shape.radius) = nan;

      P = cat(3, P1, P2, P3, P4);
      dist = cat(1, t1, t2, t3, t4);

      % Find first intersection
      [~, idx1] = min(dist, [], 1);

      % Find second intersection
      dist2 = dist;
      dist2(idx1, :) = nan;
      [~, idx2] = min(dist2, [], 1);

      % Package outputs
      dist = permute(dist(:, :, [idx1, idx2]), [1, 3, 2]);
      locs = permute(locs(:, :, [idx1, idx2]), [1, 3, 2]);
      norms = reshape(shape.normalsXyz(reshape(locs, 3, [])), 3, 2, []);
    end

    function S = surfInternal(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Converts the shape to a PatchMesh and calls surf.
      %
      % Usage
      %   S = shape.surfInternal(...)
      %
      % Optional named parameters
      %   - resolution ([ntheta, nphi]) -- Number of faces in theta/phi
      %     directions.  Default: ``[0, 20]`` uses cylinder corners for
      %     theta direction unless otherwise set.
      %
      % Additional named parameters are passed to PatchMesh.surf.

      p = inputParser;
      p.addParameter('resolution', [0, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      S = surfInternal@ott.shape.mixin.AxisymStarShape(shape, ...
          'resolution', p.Results.resolution, unmatched{:});
    end

    function shape = scaleInternal(shape, sc)
      shape.height = shape.height * sc;
      shape.radius = shape.radius * sc;
    end
  end

  methods % Getters/setters
    function r = get.maxRadius(shape)
      % Calculate the maximum particle radius
      r = sqrt(shape.radius.^2 + (shape.height/2).^2);
    end

    function v = get.volume(shape)
      % Calculate the volume
      v = pi*shape.radius.^2*shape.height;
    end

    function p = get.perimeter(shape)
      % Calculate the perimeter of the object
      p = 2.0 * (2.0*shape.radius + shape.height);
    end

    function b = get.xySymmetry(shape)
      b = true;
    end

    function bb = get.boundingBox(shape)
      bb = [[-1, 1; -1, 1].*shape.radius; [-1, 1].*shape.height./2];
    end
  end
end
