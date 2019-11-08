classdef AxisymLerp < ott.shapes.StarShape & ott.shapes.AxisymShape
%AxisymLerp a axisymmetric particle with lerping between points
% Inherits from ott.shapes.StarShape and ott.shapes.AxisymShape.
%
% properties:
%   rho        % Radial position of defining points (cylindrical coords)
%   z          % z position of defining points (cylindrical coords)
%
% See also AxisymLerp

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    rho        % Radial position of defining points (cylindrical coords)
    z          % z position of defining points (cylindrical coords)
  end

  methods
    function shape = AxisymLerp(rho, z)
      % AxisymLerp construct a rotationally symetric particle with
      % linear interpolation (lerping) between points.
      %
      % The particle profile must be star shaped.
      %
      % AxisymLerp(rho) or AxisymLerp([rho(:) z(:)]) creates
      % a new cylinder with the specified profile.

      shape = shape@ott.shapes.StarShape();

      if nargin == 2
        shape.rho = rho(:);
        shape.z = z(:);
      else
        shape.rho = rho(:, 1);
        shape.z = rho(:, 2);
      end

      % Ensure rho and z have the same size
      [shape.rho, shape.z] = ott.utils.matchsize(shape.rho, shape.z);

      % Ensure rho is positive
      shape.rho = abs(shape.rho);

      % Sort the points in ascending order from 0 to pi
      corner_angles = atan2(shape.z, shape.rho) + pi/2;
      [corner_angles, I] = sort(corner_angles);
      shape.rho = shape.rho(I);
      shape.z = shape.z(I);

      % Add extra points if the axis isn't included
      tol = shape.maxRadius * 1.0e-6;
      if abs(corner_angles(1) - 0.0) > tol
        shape.rho = [ 0; shape.rho ];
        shape.z = [ shape.z(1); shape.z ];
      end
      if abs(corner_angles(end) - pi) > tol
        shape.rho = [ shape.rho; 0 ];
        shape.z = [ shape.z; shape.z(end) ];
      end

    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = max(sqrt(shape.rho.^2 + shape.z.^2));
    end

    function v = get_volume(shape)
      % Calculate the volume of the particle
      error('Not yet implemented');
    end

    function p = get_perimiter(shape)
      % Calculate the perimiter of the object
      p = 2.0 * sum(sqrt((shape.rho(2:end) - shape.rho(1:end-1)).^2 ...
          + (shape.z(2:end) - shape.z(1:end-1)).^2));
    end

    function r = radii(shape, theta, phi)
      % RADII returns the radius for each requested point

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Corner angles should be in range [0, pi]
      corner_angles = atan2(shape.z, shape.rho) + pi/2;

      r = zeros(size(theta));

      % For each side
      for ii = 2:length(corner_angles)

        % Find points that belong to this side
        if ii == 2
          pts = theta < corner_angles(ii);
        elseif ii == length(corner_angles)
          pts = theta >= corner_angles(ii-1);
        else
          pts = theta >= corner_angles(ii-1) & theta < corner_angles(ii);
        end

        % Calculate side lengths of triangle O-C1-C2
        a = sqrt(shape.rho(ii-1).^2 + shape.z(ii-1).^2);
        b = sqrt(shape.rho(ii).^2 + shape.z(ii).^2);
        c = sqrt((shape.rho(ii) - shape.rho(ii-1)).^2 ...
            + (shape.z(ii) - shape.z(ii-1)).^2);

        % Calculate sin(angle O-C1-C2)
        kappa = b * sin(corner_angles(ii) - corner_angles(ii-1)) / c;

        % Calculate angle C1-O-P
        chi = theta(pts) - corner_angles(ii-1);

        % Calculate angle O-P-C1
        gamma = pi - chi - abs(asin(kappa));

        % Calculate and store the radius
        r(pts) = a * kappa ./ sin(gamma);
      end
    end

    function n = normals(shape, theta, phi)
      % NORMALS calculate the normals for the specified points.

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Corner angles should be in range [0, pi]
      corner_angles = atan2(shape.z, shape.rho) + pi/2;

      % Calculate angle of segments
      segment_angles = atan2(shape.z(2:end) - shape.z(1:end-1), ...
          shape.rho(2:end) - shape.rho(1:end-1)) + pi/2;

      n = zeros(size(theta, 1), 3);

      % For each side
      for ii = 2:length(corner_angles)

        % Find points that belong to this side
        if ii == 2
          pts = theta < corner_angles(ii);
        elseif ii == length(corner_angles)
          pts = theta >= corner_angles(ii-1);
        else
          pts = theta >= corner_angles(ii-1) & theta < corner_angles(ii);
        end

        % Calculate the angle between the pt and the surface normal
        psi = segment_angles(ii-1) - theta(pts) - pi/2;

        % Calculate the normal in spherical coordinates
        n(pts, 1) = cos(psi);
        n(pts, 2) = sin(psi);
      end
    end

    function varargout = axialSymmetry(shape)
      % Return the axial symmetry for the particle

      vals = [ 1, 1, 0 ];

      % If we have z-mirror symetry, x and y are 2 fold rotationally symetric
      [~, ~, zmirrorsym] = shape.mirrorSymmetry();
      if zmirrorsym
        vals(1:2) = 2;
      end

      if nargout == 1
        varargout{1} = vals;
      else
        varargout{1} = vals(1);
        varargout{2} = vals(2);
        varargout{3} = vals(3);
      end
    end

    function varargout = mirrorSymmetry(shape)
      % Return the mirror symmetry for the particle

      zmirrorsym=false;
      tol = shape.maxRadius * 1.0e-6;
      if all(abs(shape.rho - flipud(shape.rho)) < tol) ...
          && all(abs(shape.z + flipud(shape.z)) < tol)
        zmirrorsym=true;
      end

      if nargout == 1
        varargout{1} = [ true, true, zmirrorsym ];
      else
        varargout{1} = true;
        varargout{2} = true;
        varargout{3} = zmirrorsym;
      end
    end

    function varargout = surf(shape, varargin)
      % SURF generate a visualisation of the shape
      %
      % SURF() displays a visualisation of the shape in the current figure.
      %
      % [X, Y, Z] = surf() calculates the coordinates and arranges them
      % in a grid.
      %
      % Optional named arguments:
      %   npoints   num   Number of azimuthal points to use (default: 20)
      %   points    3xN   Array of points to calculate (default: [])
      %
      % See also ott.shapes.StarShape/surf

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('points', []);
      p.addParameter('npoints', 20);
      p.parse(varargin{:});

      if isempty(p.Results.points)

        theta = atan2(shape.z, shape.rho) + pi/2;
        [~, phi] = ott.utils.angulargrid(1, p.Results.npoints);

        [varargout{1:nargout}] = surf@ott.shapes.StarShape(...
            shape, varargin{:}, 'points', { theta,  phi } );
      else
        [varargout{1:nargout}] = surf@ott.shapes.StarShape(...
            shape, varargin{:} );
      end
    end
  end
end
