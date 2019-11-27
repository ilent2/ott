classdef Ellipsoid < ott.shapes.StarShape & ott.shapes.AxisymShape
%Ellipsoid a simple ellipsoid shape
%
% properties:
%   a         % x-axis scaling
%   b         % y-axis scaling
%   c         % z-axis scaling

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    a
    b
    c
  end

  methods
    function shape = Ellipsoid(a, b, c)
      % ELLIPSOID construct a new ellipsoid
      %
      % ELLIPSOID(a, b, c) or ELLIPSOID([a, b, c]) specifies the
      % scaling in the x, y and z directions.

      shape = shape@ott.shapes.StarShape();

      if nargin == 1
        shape.a = a(1);
        shape.b = a(2);
        shape.c = a(3);
      else
        shape.a = a;
        shape.b = b;
        shape.c = c;
      end

    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = max([shape.a, shape.b, shape.c]);
    end

    function v = get_volume(shape)
      % Calculate the particle volume
      v = 4./3.*pi.*shape.a.*shape.b.*shape.c;
    end

    function b = isSphere(shape)
      % ISSPHERE Returns true if the shape is a sphere
      b = shape.a == shape.b && shape.a == shape.c;
    end

    function p = get_perimiter(shape)
      % Calculate the perimiter of the object
      %
      % Only works if the shape is rotationally symmetric about z.

      % Check that the particle is axisymmetric
      axisym = shape.axialSymmetry();
      if axisym(3) ~= 0
        p = NaN;
        return;
      end

      a = max(shape.a, shape.c);
      b = min(shape.a, shape.c);

      % From mathsisfun.com/geometry/ellipse-perimiter.html
      h = (a - b).^2 ./ (a + b).^2;
      p = 1.0;
      nmax = 100;
      tol = 1.0e-3;
      nchooseks = @(n, k) gamma(n)./(gamma(n - k).*gamma(k));
      for ii = 1:nmax
        newp = nchooseks(0.5, ii).^2.*h.^ii;
        p = p + newp;
        if abs(newp) < tol
          break;
        end
      end
      p = p .* pi .* (a + b);

    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      r = 1 ./ sqrt( (cos(phi).*sin(theta)/shape.a).^2 + ...
          (sin(phi).*sin(theta)/shape.b).^2 + (cos(theta)/shape.c).^2 );
    end

    function n = normals(shape, theta, phi)
      % NORMALS calcualtes the normals for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % r = 1/sqrt( cos(phi)^2 sin(theta)^2 / a^2 +
      %     sin(phi)^2 sin(theta)^2 / b2 + cos(theta)^2 / c^2 )
      % dr/dtheta = -r^3 ( cos(phi)^2 / a^2 + sin(phi)^2 / b^2
      %             - 1/c^2 ) sin(theta) cos(theta)
      % dr/dphi = -r^3 sin(phi) cos(phi) (1/b^2 - 1/a^2) sin(theta)^2
      % sigma = rhat + thetahat r^2 sin(theta) cos(theta) *
      %     ( cos(phi)^2/a^2 + sin(phi)^2/b^2 - 1/c^2 )
      %    + phihat r^2 sin(theta) sin(phi) cos(phi) (1/b^2 - 1//a^2)

      r = shape.radii(theta, phi);
      sigma_r = ones(size(r));
      sigma_theta = r.^2 .* sin(theta) .* cos(theta) .* ...
          ( (cos(phi)/shape.a).^2 + (sin(phi)/shape.b).^2 - 1/shape.c^2 );
      sigma_phi = r.^2 .* sin(theta) .* sin(phi) .* cos(phi) .* ...
          (1/shape.b^2 - 1/shape.a^2);
      sigma_mag = sqrt( sigma_r.^2 + sigma_theta.^2 + sigma_phi.^2 );
      n = [ sigma_r./sigma_mag sigma_theta./sigma_mag ...
          sigma_phi./sigma_mag ];
    end

    function varargout = axialSymmetry(shape)
      % Return the axial symmetry for the particle

      if shape.a == shape.b; ab = 0; else; ab = 2; end
      if shape.a == shape.c; ac = 0; else; ac = 2; end
      if shape.b == shape.c; bc = 0; else; bc = 2; end

      if nargout == 1
        varargout{1} = [ bc, ac, ab ];
      else
        varargout{1} = bc;
        varargout{2} = ac;
        varargout{3} = ab;
      end
    end
  end
end
