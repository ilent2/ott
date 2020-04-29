classdef BiconcaveDisc < ott.shapes.AxisymShape
% A biconcave disc shape
%
% This shape can be used to model cells such as unstressed Red Blood Cells.
%
% Methods
%   - radialProfile -- Calculate the height as a function of radius
%
% Properties
%   - radius -- Radius of the disc
%   - a0, a1, a2 -- Coefficients describing the discs shape

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% TODO: Fix these things in version 2:
% TODO: Add a better base class
% TODO: Determine if star shaped
% TODO: Add a surf method (perhaps in base class)
% TODO: Add full list of properties and methods to header
% TODO: Add properties to constructor

  properties
    radius
    a0
    a1
    a2
  end

  methods (Hidden)
    function v = get_volume(shape)
      % Calculate volume numerically
      %
      % TODO: This should move elsewhere
      % This is common with AxisymLerp

      r = linspace(0, shape.maxRadius, 100);
      z = shape.radialProfile(r);
      v = trapz(r, z) .* 2 .* (2.*pi.*r);
    end

    function p = get_perimiter(shape)
      % Calculate perimeter numerically
      %
      % TODO: This should move elsewhere
      % This is common with AxisymLerp

      rho = linspace(0, shape.maxRadius, 100);
      z = shape.radialProfile(r);

      p = 2.0 * sum(sqrt((rho(2:end) - rho(1:end-1)).^2 ...
          + (z(2:end) - z(1:end-1)).^2));
    end

    function r = get_maxRadius(shape)
      r = shape.radius;
    end
  end

  methods
    function shape = BiconcaveDisc(radius, coefficients)
      % Construct a new biconcave disc shape representation
      %
      % Usage
      %   shape = BiconcaveDisc(radius, [a0, a1, a2])
      %
      % Parameters
      %   - radius (numeric) -- radius of the disc.
      %   - [a0, a1, a2] (numeric) -- coefficients describing the shape.
      %
      % For a shape similar to a red blood cell, use
      % ``a = [0.0518, 2.0026, -4.491]`` and ``radius = 7.82``.

      shape = shape@ott.shapes.AxisymShape();

      assert(numel(coefficients) == 3 && isnumeric(coefficients), ...
        'coefficients must be 3 element numeric array');
      assert(isnumeric(radius) && isscalar(radius), ...
        'radius must be numeric scalar');

      shape.radius = radius;
      shape.a0 = coefficients(1);
      shape.a1 = coefficients(2);
      shape.a2 = coefficients(3);
    end

    function z = radialProfile(shape, r)
      % Calculate the radial profile of the shape
      %
      % The profile for a biconcave disc is::
      %
      %   z(r) = D \sqrt{1 - \frac{4r^2}{D^2}} \left(a_0 +
      %       \frac{a_1 r^2}{D^2} + \frac{a_2 r^4}{D^4} \right)
      %
      % Usage
      %   z = shape.radialProfile(r)

      r2 = r.^2;

      z = 2.0 .* shape.radius .* sqrt(1 - r2./shape.radius.^2) ...
        .* (shape.a0 + (shape.a1.*r2) ./ shape.radius.^2 ./ 4 ...
        + (shape.a2.*r2.^2) ./ shape.radius.^4 ./ 16);
    end

    function b = inside(shape, r, t, p, varargin)

      rtp = [r(:), t(:), p(:)].';

      % Convert to xyz coordinates
      xyz = ott.utils.rtp2xyz(rtp.').';

      % Calculate surface
      rho = sqrt(sum(xyz(1:2, :).^2, 1));
      surfz = shape.radialProfile(rho);

      % Determine if point is inside
      b = abs(xyz(3, :)) < surfz & rho <= shape.maxRadius;
      b = reshape(b, size(r));
    end

    function normals(shape, theta, phi)
      % Same as lerp?
      error('not implemented yet');
    end

    function radii(shape, theta, phi)
      % Same as lerp?
      error('not implemented yet');
    end

    function varargout = surf(shape, varargin)
      % SURF generate a visualisation of the shape
      %
      % SURF() displays a visualisation of the shape in the current figure.
      %
      % [X, Y, Z] = surf() calculates the coordinates and arranges them
      % in a grid.

      r = linspace(0, shape.maxRadius, 30);
      t = linspace(0, 2*pi, 30);
      [rr, tt] = meshgrid(r, t);

      zz = shape.radialProfile(rr);

      rr = [rr, fliplr(rr)];
      tt = [tt, fliplr(tt)];
      zz = [zz, -fliplr(zz)];

      xx = rr .* sin(tt);
      yy = rr .* cos(tt);

      surf(xx, yy, zz);
      axis equal;

    end

    function s = axialSymmetry(shape)
      % Same as lerp (almost)

      s = [2, 2, 0];

      if nargout == 1
        varargout{1} = vals;
      else
        varargout{1} = vals(1);
        varargout{2} = vals(2);
        varargout{3} = vals(3);
      end
    end
  end
end

