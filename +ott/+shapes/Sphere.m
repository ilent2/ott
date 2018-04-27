classdef Sphere < ott.shapes.StarShape & ott.shapes.AxisymShape
%Sphere a simple sphere shape
%
% properties:
%   radius        % Radius of the sphere
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius
  end

  methods
    function sphere = Sphere(radius)
      % SPHERE construct a new sphere
      %
      % SPHERE(radius) specifies the radius of the sphere

      sphere = sphere@ott.shapes.StarShape();
      sphere.radius = radius;
    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = shape.radius;
    end

    function p = get_perimiter(shape)
      % Calculate the perimiter of the object
      p = 2.0 * pi * shape.radius;
    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      r = ones(size(theta)) * shape.radius;
    end

    function n = normals(shape, theta, phi)
      % NORMALS calculate the normals for the specified points.

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      n = ones(size(theta)) * [ 1, 0, 0 ];
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
end
