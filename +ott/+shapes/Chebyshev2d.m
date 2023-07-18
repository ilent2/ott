classdef Chebyshev2d < ott.shapes.StarShape & ott.shapes.AxisymShape
%Implementation of a 2-D Chebyshev particle
%
% r(theta) = radius*(1 + sscale*cos(order*theta))
%
% Properties:
%   radius        % Sphere radius
%   sscale        % Perterbation scale
%   order         % Perterbation order
%
% This implementation is based on
% github.com/michaelkahnert/Tsym-5.2

% Copyright 2023 IST Austria
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius
    sscale
    order
  end

  methods
    function shape = Chebyshev2d(radius, sscale, order)
      % CHEBYSHEV2D construct a new 2-D Chebyshev particle
      %
      % Chebyshev2d(radius, sscale, order)

      shape = shape@ott.shapes.StarShape();
      shape.radius = radius;
      shape.sscale = sscale;
      shape.order = order;
    end

    function r = get_maxRadius(shape)
      r = shape.radius + shape.sscale;
    end

    function v = get_volume(shape)
      error('Not yet implemented');
    end

    function p = get_perimiter(shape)
      error('Not yet implemented');
    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta, phi] = ott.utils.matchsize(theta, phi);

      r = shape.radius*(1.0 + shape.sscale.*cos(shape.order*theta));
    end

    function n = normals(shape, theta, phi)
      % NORMALS calculate the normals for the specified points.

      theta = theta(:);
      phi = phi(:);
      [theta, phi] = ott.utils.matchsize(theta, phi);

      rth = shape.radius*shape.order*shape.sscale*sin(shape.order*theta);
      rphi = zeros(size(theta));
      n = [ones(size(rth)), rth, rphi];
      n = n ./ vecnorm(n, 2, 2);
    end

    function varargout = axialSymmetry(shape)
      % Return the axial symmetry for the particle

      vals = [ 1, 1, 0 ];

      % If we have z-mirror symetry, x and y are 2 fold rotationally symetric
      % TODO: implement z-mirror symmetry check
      %[~, ~, zmirrorsym] = shape.mirrorSymmetry();
      %if zmirrorsym
      %  vals(1:2) = 2;
      %end

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