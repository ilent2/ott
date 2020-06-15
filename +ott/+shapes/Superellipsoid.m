classdef Superellipsoid < ott.shapes.Shape ...
    & ott.shapes.mixin.StarShape ...
    & ott.shapes.mixin.NumericalVolume ...
    & ott.shapes.mixin.IsosurfSurfPoints
% Superellipsoid shape
%
% In Cartesian coordinates, a superellipsoid is defined by::
%
%     ( |x/a|^(2/e) + |y/b|^(2/e) )^(e/n) + |z/c|^(2/n) = 1
%
% where :math:`a,b,c` are the radii along Cartesian directions,
% and :math:`e,n` are the east-west and north-south smoothness parameters.
% For more details see https://en.wikipedia.org/wiki/Superellipsoid
%
% Properties
%   - radii     -- Radii along Cartesian axes [X; Y; Z]
%   - ew        -- East-West smoothness (ew = 1 for ellipsoid)
%   - ns        -- North-South smoothness (sw = 1 for ellipsoid)

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radii      % Radii along Cartesian axes [X; Y; Z]
    ew         % East-West smoothness (ew = 1 for ellipsoid)
    ns         % North-South smoothness (sw = 1 for ellipsoid)
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    zRotSymmetry       % Degree of z rotational symmetry
    xySymmetry         % True if the particle is xy-plane mirror symmetric
  end

  methods
    function shape = Superellipsoid(varargin)
      % Construct a new superellipsoid
      %
      % Usage
      %   shape = Superellipsoid(radii, ew, ns, ...)
      %
      % Parameters
      %   - radii (3 numeric) -- Radii along Cartesian axes.
      %   - ew (numeric) -- East-west smoothness (xy-plane)  Default: 0.8
      %   - ns (numeric) -- North-south smoothness (z-axis)  Default: 1.2

      p = inputParser;
      p.addOptional('radii', [1, 2, 3]);
      p.addOptional('ew', 0.8);
      p.addOptional('ns', 1.2);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shapes.Shape(unmatched{:});
      shape.radii = p.Results.radii;
      shape.ew = p.Results.ew;
      shape.ns = p.Results.ns;
    end

    function r = starRadii(shape, theta, phi)
      % Return the ellipsoid radii at specified angles
      %
      % Usage
      %   r = shape.starRadii(theta, phi)

      assert(all(size(theta) == size(phi)), ...
          'theta and phi must have same size');

      acp = abs(cos(phi));
      asp = abs(sin(phi));
      act = abs(cos(theta));
      ast = abs(sin(theta));
      cpsp = (acp/shape.radii(1)).^(2/shape.ew) ...
          + (asp/shape.radii(2)).^(2/shape.ew);
      r = ( ast.^(2/shape.ns) .* cpsp.^(shape.ew/shape.ns) ...
          + (act/shape.radii(3)).^(2/shape.ns) ).^(-shape.ns/2);
    end

    function n = normalsRtpInternal(shape, rtp)
      % Calculate the normals for each requested point

      theta = rtp(2, :).';
      phi = rtp(3, :).';

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

      n = [ sigma_r./sigma_mag, sigma_theta./sigma_mag, ...
         sigma_phi./sigma_mag ].';

      n = ott.utils.rtpv2xyzv(n, rtp);
    end
  end

  methods % Getters/setters
    function shape = set.radii(shape, val)
      assert(isnumeric(val) && numel(val) == 3 && all(val >= 0), ...
          'radii must be positive numeric 3-vector');
      shape.radii = val(:);
    end

    function r = get.maxRadius(shape)
      % Calculate the maximum particle radius

      % This is likely an over-estimate for most particles, the lower
      % bounds of this would be max(shape.radii).
      r = vecnorm(shape.radii);
    end

    function bb = get.boundingBox(shape)
      bb = [-1, 1; -1, 1; -1, 1].*shape.radii;
    end

    function b = get.xySymmetry(shape)
      b = true;
    end
    function q = get.zRotSymmetry(shape)
      if shape.radii(1) == shape.radii(2)
        if shape.ew == 1.0
          q = 0;
        else
          q = 4;
        end
      else
        q = 2;
      end
    end
  end
end
