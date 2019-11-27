classdef RectangularPrism < ott.shapes.StarShape
%Cube a simple cube shape
%
% properties:
%   x        Size of prism in x direction
%   y        Size of prism in y direction
%   z        Size of prism in z direction
%
% See also RectangularPrism and ott.shapes.Cube.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    x     % Size of prism in x direction
    y     % Size of prism in y direction
    z     % Size of prism in z direction
  end

  methods
    function shape = RectangularPrism(x, y, z)
      % RectangularPrism construct a new Rectangular Prism
      %
      % RectangularPrism(x, y, z) construct a new rectangular prism
      % with edges lengths specified by x, y, and z.

      shape = shape@ott.shapes.StarShape();
      shape.x = x;
      shape.y = y;
      shape.z = z;
    end

    function r = get_maxRadius(shape)
      % Calculate the maximum particle radius
      r = sqrt(sum([shape.x, shape.y, shape.z]/2.^2));
    end

    function v = get_volume(shape)
      % Calculate the volume of the cube
      v = shape.x*shape.y*shape.z;
    end

    function [f, xyz] = faces(shape, theta, phi)
      % FACE determine which face a point is on
      % Edges and corners belong to multiple faces
      
      % Coordinate transform (squished coords)
      theta = acos(shape.z*cos(theta)./sqrt(...
        shape.x^2.*cos(phi).^2.*sin(theta).^2 + ...
        shape.y^2.*sin(phi).^2.*sin(theta).^2 + ...
        shape.z^2.*cos(theta).^2));
      phi = atan(tan(phi)*shape.x/shape.y);

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

    function varargout = locations(shape, theta, phi)
      % LOCATIONS calculate the Cartessian coordinate locations

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Get the face this point belongs to
      [planes, xyz] = shape.faces(theta, phi);
      xplane = logical(planes(:, 1));
      yplane = logical(planes(:, 2));
      zplane = logical(planes(:, 3));

      % For points, ensure the point belongs to only one plane
      yplane = yplane & ~xplane;
      zplane = zplane & ~xplane & ~yplane;

      % Rescale the vectors so they sit on the cube
      xyz(xplane, :) = xyz(xplane, :) ./ abs(xyz(xplane, 1));
      xyz(yplane, :) = xyz(yplane, :) ./ abs(xyz(yplane, 2));
      xyz(zplane, :) = xyz(zplane, :) ./ abs(xyz(zplane, 3));

      xyz = xyz .* [shape.x, shape.y, shape.z] / 2.0;

      if nargout == 1
        varargout{1} = xyz;
      else
        varargout{1} = xyz(:, 1);
        varargout{2} = xyz(:, 2);
        varargout{3} = xyz(:, 3);
      end
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
      
      sym = [ shape.y == shape.z, shape.x == shape.z, shape.x == shape.y];

      if nargout == 1
        varargout{1} = sym*4;
      else
        varargout{1} = sym(1)*4;
        varargout{2} = sym(2)*4;
        varargout{3} = sym(3)*4;
      end
    end

    function varargout = surf(shape, varargin)
      % SURF generate a visualisation of the shape
      %
      % SURF(...) displays a visualisation of the shape in the current figure.
      %
      % [X, Y, Z] = surf() calculates the coordinates and arranges them
      % in a grid suitable for use with matlab surf function.
      %
      % Optional named arguments:
      %   offset   [x;y;z]   offset for location of surface
      %   rotation   mat     rotation matrix to apply to surface
      %   points   { theta, phi }  specify points to use for surface
      %   surfoptions   {varargin} options to be passed to surf.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('points', []);
      p.addParameter('npoints', [100, 100]);
      p.addParameter('noendcap', false);
      p.addParameter('surfoptions', {});
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('axes', []);
      p.parse(varargin{:});

      if isempty(p.Results.points)
        
       verts = [ 1, 1, 1; -1, 1, 1; ...
                1, -1, 1; -1, -1, 1; ...
                1, 1, -1; -1, 1, -1; ...
                1, -1, -1; -1, -1, -1].' .* [shape.x; shape.y; shape.z]./2;
              
        idx = [1, 2, 4, 3, 1];
        idx = [idx; idx; idx+4; idx+4];
        
        X = reshape(verts(1, idx), size(idx));
        Y = reshape(verts(2, idx), size(idx));
        Z = reshape(verts(3, idx), size(idx));
        
        % Add end caps
        if p.Results.noendcap
          X(1, :) = []; X(3, :) = [];
          Y(1, :) = []; Y(3, :) = [];
        else
          X(1, :) = 0; X(4, :) = 0;
          Y(1, :) = 0; Y(4, :) = 0;
        end
        
        % Apply a rotation to the surface
        if ~isempty(p.Results.rotation)
          XYZ = [X; Y; Z];
          XYZ = p.Results.rotation * XYZ;
          X = XYZ(1, :);
          Y = XYZ(2, :);
          Z = XYZ(3, :);
        end

        % Apply a translation to the surface
        if ~isempty(p.Results.position)
          X = X + p.Results.position(1);
          Y = Y + p.Results.position(2);
          Z = Z + p.Results.position(3);
        end
        
        % Generate the surface
        if nargout == 0 || ~isempty(p.Results.axes)

          % Place the surface in the specified axes
          our_axes = p.Results.axes;
          if isempty(our_axes)
            our_axes = axes();
          end

          surf(our_axes, X, Y, Z, p.Results.surfoptions{:});
        end

        % Set outputs if requested
        if nargout ~= 0
          varargout = { X, Y, Z };
        end
      else
        [varargout{1:nargout}] = surf@ott.shapes.StarShape(...
            shape, varargin{:});
      end
    end

    function writeWavefrontObj(shape, filename)
      % Write representation of shape to Wavefront OBJ file
      %
      % writeWavefrontObj(filename) writes the shape to the given file.

      % Generate array of vertices
      verts = [ 1, 1, 1; 1, 1, -1; ...
                1, -1, 1; 1, -1, -1; ...
                -1, 1, 1; -1, 1, -1; ...
                -1, -1, 1; -1, -1, -1].' .* [shape.x; shape.y; shape.z]./2;

      % Generate array of faces
      % Order vertices so normals face outwards
      faces = [ 1, 5, 7, 3; 2, 4, 8, 6; ...
                1, 3, 4, 2; 3, 7, 8, 4; ...
                7, 5, 6, 8; 5, 1, 2, 6 ].';

      % Write the file
      shape.writeWavefrontObj_helper(filename, verts, faces);
    end
  end
end
