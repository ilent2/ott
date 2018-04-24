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

    function varargout = mirrorSymmetry(shape)
      % Return the mirror symmetry for the particle

      if nargout == 1
        varargout{1} = [ true, true, true ];
      else
        varargout{1} = true;
        varargout{2} = true;
        varargout{3} = true;
      end
    end

    function varargout = angulargrid(shape, Nmax)
      % ANGULARGRID calculate the angular grid and radii for the shape
      %
      % This is the default function with no symmetry optimisations.
      %
      % [theta, phi] = ANGULARGRID(Nmax) gets the default angular
      % grid for the particle.
      %
      % rtp = ANGULARGRID(Nmax) or [r, theta, phi] = ANGULARGRID(Nmax)
      % calculate the radii for the locations theta, phi.

      ntheta = 4*(Nmax + 2);
      nphi = 1;

      % TODO: There is also mirror symmetry

      % Calculate the angular grid
      [theta, phi] = ott.utils.angulargrid(ntheta, nphi);

      if nargout == 2
        varargout{1} = theta;
        varargout{2} = phi;
      else
        % Calculate the radii
        r = shape.radii(theta, phi);
        if nargout == 1
          varargout{1} = [ r theta phi ];
        else
          varargtou{1} = r;
          varargtou{2} = theta;
          varargtou{3} = phi;
        end
      end
    end
  end
end
