classdef Cube < ott.shapes.StarShape
%Cube a simple cube shape
%
% properties:
%   width        % Width of the cube
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    width
  end

  methods
    function shape = Cube(width)
      % CUBE construct a new cube
      %
      % CUBE(width) specifies the width of the cube

      shape = shape@ott.shapes.StarShape();
      shape.width = width;
    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = sqrt(3*(shape.width/2).^2);
    end

    function [f, xyz] = faces(shape, theta, phi)
      % FACE determine which face a point is on
      % Edges and corners belong to multiple faces

      % First find the rough Cartesian coordinates
      xyz = ott.utils.rtp2xyz(sqrt(3), theta, phi);

      % Round coordinates towards zero (should all now be 0, 1 or -1)
      xyz_sign = fix(xyz);

      % Determine which face the point is on
      xloop = xyz_sign(:, 1) == 0;
      yloop = xyz_sign(:, 2) == 0;
      zloop = xyz_sign(:, 3) == 0;

      xplane = ~xloop & abs(xyz(:, 1)) >= abs(xyz(:, 2)) ...
          & abs(xyz(:, 1)) >= abs(xyz(:, 3));
      yplane = ~yloop & abs(xyz(:, 2)) >= abs(xyz(:, 1)) ...
          & abs(xyz(:, 2)) >= abs(xyz(:, 3));
      zplane = ~zloop & abs(xyz(:, 3)) >= abs(xyz(:, 1)) ...
          & abs(xyz(:, 3)) >= abs(xyz(:, 2));
      corners = ~xplane & ~yplane & ~zplane;

      % Store the result
      f = zeros(size(xyz_sign));
      f(xplane, 1) = xyz_sign(xplane, 1);
      f(yplane, 2) = xyz_sign(yplane, 2);
      f(zplane, 3) = xyz_sign(zplane, 3);
      f(corners, :) = sign(xyz(corners, :));
    end

    function nxyz = normalsXyz(shape, theta, phi)
      % NORMALSXYZ calculates Cartessian normals

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Determine which face we are on (i.e. the Cartesian normals)
      nxyz = shape.faces(theta, phi);

      % Normalize the normals
      nxyz = nxyz ./ sqrt(dot(nxyz, nxyz, 2));
    end

    function xyz = locations(shape, theta, phi)
      % LOCATIONS calculate the Cartessian coordinate locations

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Get the face this point belongs to
      [planes, xyz] = shape.faces(theta, phi);
      xplane = logical(planes(:, 1));
      yplane = logical(planes(:, 2));
      zplane = logical(planes(:, 3));

      % Rescale the vectors so they sit on the cube
      xyz(xplane, :) = xyz(xplane, :) ./ abs(xyz(xplane, 1));
      xyz(yplane, :) = xyz(yplane, :) ./ abs(xyz(yplane, 2));
      xyz(zplane, :) = xyz(zplane, :) ./ abs(xyz(zplane, 3));

      xyz = shape.width .* xyz / 2.0;
    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point
      [r, ~, ~] = ott.utils.xyz2rtp(shape.locations(theta, phi));
    end

    function n = normals(shape, theta, phi)
      % NORMALS calculate the normals for the specified points.

      % Calculate the Cartessian normals
      nxyz = shape.normalsXyz(theta, phi);

      % Convert from cartesian to spherical normals
      n = ott.utils.xyzv2rtpv(nxyz, shape.locations(theta, phi));
    end

    function varargout = axialSymmetry(shape)
      % Return the axial symmetry for the particle

      if nargout == 1
        varargout{1} = [ 4, 4, 4 ];
      else
        varargout{1} = 4;
        varargout{2} = 4;
        varargout{3} = 4;
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
  end
end
