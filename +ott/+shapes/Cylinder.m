classdef Cylinder < ott.shapes.StarShape & ott.shapes.AxisymShape
%Cylinder a simple cylinder shape
%
% properties:
%   radius        % Radius of the cylinder
%   height        % Height of the cylinder
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius        % Radius of the cylinder
    height        % Height of the cylinder
  end

  methods
    function shape = Cylinder(radius, height)
      % CYLINDER construct a new cylinder
      %
      % CYLINDER(radius, height) or CYLINDER([radius, height]) creates
      % a new cylinder with the specified radius and height.

      shape = shape@ott.shapes.StarShape();

      if nargin == 2
        shape.radius = radius;
        shape.height = height;
      else
        shape.radius = radius(1);
        shape.height = radius(2);
      end
    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = sqrt(shape.radius.^2 + (shape.height/2).^2);
    end

    function p = get_perimiter(shape)
      % Calculate the perimiter of the object
      p = 2.0 * (2.0*shape.radius + shape.height);
    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

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

    function n = normals(shape, theta, phi)
      % NORMALS calculate the normals for the specified points.

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Find angle of edges from xy plane
      edge_angle = atan(shape.height/(2*shape.radius));

      sides = ( theta >= pi/2 - edge_angle ) ...
          & ( theta <= pi/2 + edge_angle ) ;
      top = theta < pi/2 - edge_angle;
      bottom = theta > pi/2 + edge_angle;
      ends = top | bottom;

      n = zeros(size(theta, 1), 3);
      n(sides,1) = sin(theta(sides));
      n(sides,2) = cos(theta(sides));
      n(ends,1) = abs(cos(theta(ends)));
      n(ends,2) = - sin(theta(ends)) .* sign(cos(theta(ends)));
    end

    function varargout = axialSymmetry(shape)
      % Return the axial symmetry for the particle

      if nargout == 1
        varargout{1} = [ 2, 2, 0 ];
      else
        varargout{1} = 2;
        varargout{2} = 2;
        varargout{3} = 0;
      end
    end
  end
end
