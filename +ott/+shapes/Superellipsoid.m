classdef Superellipsoid < ott.shapes.StarShape
%Superellipsoid a simple superellipsoid shape
%
% properties:
%   a         % x-axis scaling
%   b         % y-axis scaling
%   c         % z-axis scaling
%   ew        % East-West smoothness (ew = 1 for ellipsoid)
%   ns        % North-South smoothness (sw = 1 for ellipsoid)

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    a
    b
    c
    ew
    ns
  end

  methods
    function shape = Superellipsoid(a, b, c, ew, ns)
      % SUPERELLIPSOID construct a new superellipsoid
      %
      % SUPERELLIPSOID(a, b, c, e, n) or SUPERELLIPSOID([a, b, c, e, n])
      % specifies the scaling in the x, y and z directions.
      %
      % Defined by (cartesian)
      %   { (x/a)^(2/e) + (y/b)^(2/e) }^(e/n) + (z/c)^(2/n) = 1
      % where n = NS roundedness, e = EW roundedness

      shape = shape@ott.shapes.StarShape();

      if nargin == 1
        shape.a = a(1);
        shape.b = a(2);
        shape.c = a(3);
        shape.ew = a(4);
        shape.ns = a(5);
      else
        shape.a = a;
        shape.b = b;
        shape.c = c;
        shape.ew = ew;
        shape.ns = ns;
      end

    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = max([shape.a, shape.b, shape.c]);
    end

    function v = get_volume(shape)
      % Calculate the volume of the particle
      error('Not yet implemented');
    end

    function b = isSphere(shape)
      % ISSPHERE Returns true if the shape is a sphere
      b = shape.a == shape.b && shape.a == shape.c ...
          && shape.ew == 1 && shape.ns == 1;
    end

    function b = isEllipsoid(shape)
      % ISSPHERE Returns true if the shape is a sphere
      b = shape.ew == 1 && shape.ns == 1;
    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      acp = abs(cos(phi));
      asp = abs(sin(phi));
      act = abs(cos(theta));
      ast = abs(sin(theta));
      cpsp = (acp/shape.a).^(2/shape.ew) + (asp/shape.b).^(2/shape.ew);
      r = ( ast.^(2/shape.ns) .* cpsp.^(shape.ew/shape.ns) ...
          + (act/shape.c).^(2/shape.ns) ).^(-shape.ns/2);
    end

    function n = normals(shape, theta, phi)
      % NORMALS calcualtes the normals for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Converting to spherical coordinates,
      % r = [ sin(theta)^(2/n) { (cos(phi)/a)^(2/e)
      %     + (sin(phi)/b)^(2/e) }^(e/n)
      %     + (cos(theta)/c)^(2/n) ]^(-n/2)
      % dr/dtheta = -r^((2+n)/n) sin(theta) cos(theta)
      %   x [ sin(theta)^((1-n)/n) { (cos(phi)/a)^(2/e)
      %             + (sin(phi)/b)^(2/e) }^(e/n)
      %             - cos(theta)^((1-n)/n) * 1/c^(2/n) )
      % dr/dphi = -r^((2+n)/n) sin(phi) cos(phi)
      %   x { (cos(phi)/a)^(2/e) + (sin(phi)/b)^(2/e) }^((e-n)/n)
      %   x ( sin(phi)^((1-e)/e)/b^(2/e)
      %         - cos(phi)^((1-e)/e)/a^(2/e)) sin(theta)^(2/n)
      % sigma(theta,phi) = rhat - 1/r dr/dtheta * thetahat
      %                    - 1/(r sin(theta)) * dr/dphi * phihat

      acp = abs(cos(phi));
      asp = abs(sin(phi));
      act = abs(cos(theta));
      ast = abs(sin(theta));
      cpsp = (acp/shape.a).^(2/shape.ew) + (asp/shape.b).^(2/shape.ew);

      r = shape.radii(theta, phi);

      sigma_r = ones(size(r));
      sigma_theta = r.^(2/shape.ns) .* ast .* cos(theta) .* ...
          ( ast.^((1-shape.ns)/shape.ns) .* cpsp.^(shape.ew/shape.ns) ...
          - act.^((1-shape.ns)/shape.ns)/shape.c^(2/shape.ns) );
      sigma_phi = r.^(2/shape.ns) .* ast.^((2-shape.ns)/shape.ns) ...
          .* sin(phi) .* cos(phi) .* cpsp.^((shape.ew-shape.ns)/shape.ns) ...
          .* ( asp.^((1-shape.ew)/shape.ew)/shape.b^(2/shape.ew) ...
          - acp.^((1-shape.ew)/shape.ew)/shape.a^(2/shape.ew) );
      sigma_mag = sqrt( sigma_r.^2 + sigma_theta.^2 + sigma_phi.^2 );

      n = [ sigma_r./sigma_mag sigma_theta./sigma_mag ...
         sigma_phi./sigma_mag ];
    end

    function varargout = axialSymmetry(shape)
      % Return the axial symmetry for the particle

      % TODO: Does this shape need a mirrorSymmetry function?

      % TODO: Fix these up
      bc = 1;
      ac = 1;
      ab = 1;
      if shape.ew == 1 && shape.a == shape.b
        ab = 0;
        if shape.ns == 1 && shape.a == shape.c
          bc = 0;
          ac = 0;
        end
      end

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
