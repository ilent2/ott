classdef Sphere < ott.shapes.StarShape & ott.shapes.AxisymShape
% Spherical particle shape.
% Inherits from :class:`StartShape` and :class:`AxisymShape`.
%
% properties:
%   radius        % Radius of the sphere

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius
  end

  methods
    function sphere = Sphere(radius, varargin)
      % SPHERE construct a new sphere
      %
      % Usage
      %   shape = Sphere(radius, ...)
      %
      % Optional named arguments
      %   - position (3xN numeric) -- Position of the plane.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3N numeric) -- Plane orientations.
      %     Default: ``eye(3)``.

      sphere = sphere@ott.shapes.StarShape(varargin{:});
      sphere = sphere@ott.shapes.AxisymShape();

      Nradius = numel(radius);
      assert(numel(sphere) == 1 || Nradius == 1 ...
          || Nradius == numel(sphere), ...
          'length of radius must match length of position');

      if numel(sphere) == 1 && Nradius ~= 1
        sphere = repelem(sphere, 1, Nradius);
      end

      % Assign radius
      if Nradius > 1
        for ii = 1:Nradius
          sphere(ii).radius = radius(ii);
        end
      else
        [sphere.radius] = deal(radius);
      end
    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

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

    function varargout = axialSymmetry(shape)
      % Return the axial symmetry for the particle

      if nargout == 1
        varargout{1} = [ 0, 0, 0 ];
      else
        varargout{1} = 0;
        varargout{2} = 0;
        varargout{3} = 0;
      end
    end
  end

  methods (Hidden)
    function nxyz = normalsRtpInternal(shape, rtp)
      % Calculate normals at the specified surface locations

      nxyz = ott.utils.rtpv2xyzv(repmat([1;0;0], 1, size(rtp, 2)), rtp);
    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = shape.radius;
    end

    function v = get_volume(shape)
      % Calculate the volume
      v = 4./3.*pi.*shape.radius.^3;
    end

    function p = get_perimiter(shape)
      % Calculate the perimiter of the object
      p = 2.0 * pi * shape.radius;
    end
  end

  methods % Setters/getters
    function shape = set.radius(shape, val)
      assert(isnumeric(val) && isscalar(val), ...
          'radius must be numeric scalar');
      shape.radius = val;
    end
  end
end
