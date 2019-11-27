classdef Cylinder < ott.shapes.StarShape & ott.shapes.AxisymShape
%Cylinder a simple cylinder shape
%
% properties:
%   radius        % Radius of the cylinder
%   height        % Height of the cylinder

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

    function v = get_volume(shape)
      % Calculate the volume
      v = pi*shape.radius.^2*shape.height;
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

    function [rtp, n, ds] = boundarypoints(shape, varargin)
      % BOUNDARYPOINTS calculates boundary points for surface integral
      %
      % [rtp, n, ds] = BOUDNARYPOINTS(npts) calculates the boundary points
      % and surface normal vectors in spherical coordinates and the area
      % elements of each ring.
      %
      % BOUNDARYPOINTS('Nmax', Nmax) takes a guess at a suitable npts
      % for the given Nmax.

      rho = [0.0; 1.0; 1.0; 0.0].*shape.radius;
      z = [-0.5; -0.5; 0.5; 0.5].*shape.height;

      [rtp, n, ds] = shape.boundarypoints_rhoz(rho, z, varargin{:});
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

    function varargout = surf(shape, varargin)
      % SURF generate a visualisation of the shape
      %
      % SURF() displays a visualisation of the shape in the current figure.
      %
      % [X, Y, Z] = surf() calculates the coordinates and arranges them
      % in a grid.

      p = inputParser;
      
      % Options for StarShape
      p.addParameter('points', []);
      p.addParameter('npoints', 20);    % changed from StarShape
      p.addParameter('surfoptions', {});
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('axes', []);
      
      p.addParameter('noendcap', false);
      p.parse(varargin{:});

      if isempty(p.Results.points)

        % Calculate the locations of the interesting points
        edge_angle = atan(shape.height/(2*shape.radius));
        theta = [ 0, pi/2 - edge_angle, pi/2 + edge_angle, pi ];
        [~, phi] = ott.utils.angulargrid(1, p.Results.npoints(1));

        % Remove endcaps if requested
        if p.Results.noendcap
          theta = theta(2:end-1);
        end

        [varargout{1:nargout}] = surf@ott.shapes.StarShape(shape, ...
            'points', { theta,  phi }, ...
            'npoints', [], ...
            'surfoptions', p.Results.surfoptions, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, ...
            'axes', p.Results.axes);
      else
        [varargout{1:nargout}] = surf@ott.shapes.StarShape(shape, ...
            'point', p.Results.points, ...
            'npoints', [], ...
            'surfoptions', p.Results.surfoptions, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, ...
            'axes', p.Results.axes);
      end
    end
  end
end
