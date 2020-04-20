classdef StarShape < ott.shapes.Shape
%StarShape abstract class for star shaped particles
%
% Abstract methods:
%   radii           Calculates the particle radii for angular coordinates
%   normals         Calculates the particle normals for angular coorindates
%   axialSymmetry   Returns x, y, z rotational symmetry (0 for infinite)
%   mirrorSymmetry  Returns x, y, z mirror symmetry

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Abstract)
    radii(shape, theta, phi);
    normals(shape, theta, phi);
    axialSymmetry(shape);
  end

  methods

    function varargout = mirrorSymmetry(shape)
      % Return the mirror symmetry for the particle
      %
      % Tries to determine the objects mirror symmetry from the
      % axialSymmetry.  Should be overridden for more complex objects.

      % First calculate the axial symmetry
      axialSym = shape.axialSymmetry();
      orthSym = mod(axialSym, 2) == 0;
      mirrorSym = [ orthSym(2) | orthSym(3), ...
          orthSym(1) | orthSym(3), orthSym(1) | orthSym(2) ];

      if nargout == 1
        varargout{1} = mirrorSym;
      else
        varargout{1} = mirrorSym(1);
        varargout{2} = mirrorSym(2);
        varargout{3} = mirrorSym(3);
      end
    end

    function varargout = locations(shape, theta, phi)
      % LOCATIONS calculates Cartessian coordinate locations for points
      %
      % [x, y, z] = locations(theta, phi) calculates the Cartesian coordinates
      % for locations on the surface specified by polar angle theta [0, pi)
      % and azimuthal angle phi [0, 2*pi).
      %
      % xyz = locations(theta, phi) as above, but returns the output
      % into a Nx3 matrix [x, y, z].

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      [varargout{1:nargout}] = ott.utils.rtp2xyz(...
          shape.radii(theta, phi), theta, phi);
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
      %   position   [x;y;z]   offset for location of surface
      %   rotation   mat     rotation matrix to apply to surface
      %   points   { theta, phi }  specify points to use for surface
      %   npoints  [ntheta, nphi]  specify number of points in each direction
      %   axes       []        axis to place surface in (default: gca)
      %   surfoptions   {varargin} options to be passed to surf.

      p = inputParser;
      p.addParameter('points', []);
      p.addParameter('npoints', [100, 100]);
      p.addParameter('surfoptions', {});
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('axes', []);
      p.parse(varargin{:});

      % Get the points to use for the surface
      if isempty(p.Results.points)

        % Get the size from the user inputs
        sz = p.Results.npoints;
        if numel(sz) == 1
          sz = [sz sz];
        end

        [theta, phi] = shape.angulargrid('full', true, 'size', sz);
      else
        theta = p.Results.points{1};
        phi = p.Results.points{2};

        if min(size(theta)) == 1 && min(size(phi)) == 1
          [phi, theta] = meshgrid(phi, theta);
        elseif size(theta) ~= size(phi)
          error('theta and phi must be vectors or matricies of the same size');
        end

        sz = size(theta);
      end

      % Calculate Cartesian coordinates
      [X, Y, Z] = shape.locations(theta, phi);

      % Apply a rotation to the surface
      if ~isempty(p.Results.rotation)
        XYZ = [X, Y, Z].';
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

      % Reshape the output to match the size
      X = reshape(X, sz);
      Y = reshape(Y, sz);
      Z = reshape(Z, sz);

      % Complete the sphere (add the missing faces)
      X(:, end+1) = X(:, 1);
      Y(:, end+1) = Y(:, 1);
      Z(:, end+1) = Z(:, 1);

      % Generate the surface
      if nargout == 0 || ~isempty(p.Results.axes)
        
        % Place the surface in the specified axes
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = gca();
        end
        
        surf(our_axes, X, Y, Z, p.Results.surfoptions{:});
      end
      
      % Set outputs if requested
      if nargout ~= 0
        varargout = { X, Y, Z };
      end
    end

    function varargout = voxels(shape, spacing, varargin)
      % Generate an array of xyz coordinates for voxels inside the shape
      %
      % voxels(spacing) shows a visualisation of the shape with
      % circles placed at locations on a Cartesian grid.
      %
      % xyz = voxels(spacing) returns the voxel locations.
      %
      % Optional named arguments:
      %   - 'plotoptions'   Options to pass to the plot3 function
      %   - 'visualise'     Show the visualisation (default: nargout == 0)
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.
      %     Default: 'shape'.
      %   - even_range (logical) -- Ensure the number of dipoles along
      %     any axis is even, so to avoid 0.  For example:
      %     When true, a range could be ``[-1.5, -0.5, 0.5, 1.5]``.
      %     When false, a range might be ``[-1, 0, 1]``.
      %     Default: ``false``.

      p = inputParser;
      p.addParameter('plotoptions', {...
          'MarkerFaceColor', 'w', ...
          'MarkerEdgeColor', [.5 .5 .5], ...
          'MarkerSize', 20*spacing/shape.maxRadius});
      p.addParameter('visualise', nargout == 0);
      p.addParameter('origin', 'shape');
      p.addParameter('even_range', false);
      p.parse(varargin{:});

      % Calculate range of dipoles
      numr = ceil(shape.maxRadius / spacing);
      
      if p.Results.even_range
        % Add an extra point so we don't have a point around zero
        numr = numr + 0.5;
      end
      rrange = (-numr:numr)*spacing;

      % Generate the voxel grid
      [xx, yy, zz] = meshgrid(rrange, rrange, rrange);

      % Determine which points are inside
      mask = shape.insideXyz(xx, yy, zz, 'origin', 'shape');
      xyz = [xx(mask).'; yy(mask).'; zz(mask).'];

      % Translate to world origin
      if strcmpi(p.Results.origin, 'world')
        xyz = xyz + shape.position;
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end

      % Visualise the result
      if p.Results.visualise
        plot3(xyz(1,:), xyz(2,:), xyz(3,:), 'o', p.Results.plotoptions{:});
        axis equal
        title(['spacing = ' num2str(spacing) ', N = ' int2str(sum(mask))])
      end

      % Assign output
      if nargout ~= 0
        varargout = {xyz};
      end
    end

    function n = normalsXyz(shape, theta, phi)
      % NORMALSXYZ calculates Cartessian normals
      %
      % n = normalsXyz(theta, phi) calculates the normals to the
      % surface in Cartesian coordinates and returns a Nx3 matrix.

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Convert the normals to cartesian
      n = ott.utils.rtpv2xyzv(shape.normals(theta, phi), ...
          [ zeros(size(theta)), theta, phi ]);
    end

    function b = inside(shape, radius, theta, phi, varargin)
      % INSIDE determine if point is inside the shape
      %
      % b = inside(shape, radius, theta, phi) determine if the
      % point described by radius, theta (polar), phi (azimuthal)
      % is inside the shape.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.

      p = inputParser;
      p.addParameter('origin', 'world');
      p.parse(varargin{:});

      theta = theta(:);
      phi = phi(:);
      radius = radius(:);
      [radius,theta,phi] = ott.utils.matchsize(radius,theta,phi);

      % Translate to shape origin
      if strcmpi(p.Results.origin, 'world')

        % Only do work if we need to
        if vecnorm(shape.position) ~= 0
          [x,y,z] = ott.utils.rtp2xyz(radius, theta, phi);
          x = x - shape.position(1);
          y = y - shape.position(2);
          z = z - shape.position(3);
          [radius, theta, phi] = ott.utils.xyz2rtp(x, y, z);
        end
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end

      assert(all(radius >= 0), 'Radii must be positive');

      % Determine if points are less than shape radii
      r = shape.radii(theta, phi);
      b = radius < r;

    end

    function b = insideXyz(shape, varargin)
      % INSIDEXYZ determine if Cartesian point is inside the shape
      %
      % Usage
      %   b = inside(shape, x, y, z) determine if the Cartesian point
      %   [x, y, z] is inside the star shaped object.
      %
      %   b = insideXyz(shape, xyz) as above, but using a 3xN matrix of
      %   [x; y; z] positions.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.
      %
      % See also INSIDE.

      p = inputParser;
      p.addOptional('x', []);
      p.addOptional('y', []);
      p.addOptional('z', []);
      p.addParameter('origin', 'world');
      p.parse(varargin{:});
      
      % TODO: This should be a function and added to other shapes
      if isempty(p.Results.y) && isempty(p.Results.z)
        x = p.Results.x(1, :);
        y = p.Results.x(2, :);
        z = p.Results.x(3, :);
      else
        x = p.Results.x(:);
        y = p.Results.y(:);
        z = p.Results.z(:);
        [x, y, z] = ott.utils.matchsize(x, y, z);
      end

      % Translate to shape origin
      % TODO: This should be a function and added to other shapes
      if strcmpi(p.Results.origin, 'world')
        x = x - shape.position(1);
        y = y - shape.position(2);
        z = z - shape.position(3);
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end

      % Convert to spherical coordinates
      [r, t, p] = ott.utils.xyz2rtp(x, y, z);

      % Call the spherical coordinate version
      b = shape.inside(r, t, p, 'origin', 'shape');
    end

    function varargout = angulargrid(shape, varargin)
      % ANGULARGRID calculate the angular grid and radii for the shape
      %
      % This is the default function with no symmetry optimisations.
      %
      % [theta, phi] = ANGULARGRID(Nmax) gets the default angular
      % grid for the particle.
      %
      % rtp = ANGULARGRID(Nmax) or [r, theta, phi] = ANGULARGRID(Nmax)
      % calculate the radii for the locations theta, phi.
      %
      % ANGULARGRID() uses a default Nmax of 100.
      %
      % ANGULARGRID(..., 'full', full) calculates
      % an angular grid over the full sphere.
      %
      % ANGULARGRID(..., 'size', [ntheta, nphi]) uses ntheta and
      % nphi instead of Nmax for angular grid.

      p = inputParser;
      p.addOptional('Nmax', 100);
      p.addParameter('full', false);    % Not used
      p.addParameter('size', []);    % Not used
      p.parse(varargin{:});

      if isempty(p.Results.size)

        % OTTv1 used something like the following, use it for now
        % until we can think of something better.

        ntheta = 2*(p.Results.Nmax + 2);
        nphi = 3*(p.Results.Nmax + 2) + 1;

        if ~p.Results.full

          [~, ~, z_axial_symmetry] = shape.axialSymmetry();
          if z_axial_symmetry == 0
            ntheta = 4*(p.Results.Nmax + 2);
            nphi = 1;
          else
            nphi = round(nphi / z_axial_symmetry);
          end

          [~, ~, z_mirror_symmetry] = shape.mirrorSymmetry();
          if z_mirror_symmetry
            ntheta = round(ntheta / 2);
          end
        end
      else
        ntheta = p.Results.size(1);
        nphi = p.Results.size(2);
      end

      % Special case for inifite axial symmetry
      [~, ~, z_axial_symmetry] = shape.axialSymmetry();
      if ~p.Results.full && z_axial_symmetry == 0
        nphi = 1;
      end

      % Calculate the angular grid
      [theta, phi] = ott.utils.angulargrid(ntheta, nphi);

      % Reduce the grid using z-symmetry and mirror symmetry
      if ~p.Results.full
        [~, ~, z_mirror_symmetry] = shape.mirrorSymmetry();
        if z_mirror_symmetry
          theta = theta / 2.0;    % [0, pi] -> [0, pi/2]
        end

        if z_axial_symmetry > 1
          phi = phi / z_axial_symmetry;  % [0, 2pi] -> [0, 2pi/p]
        end
      end

      if nargout == 2
        varargout{1} = theta;
        varargout{2} = phi;
      else
        % Calculate the radii
        r = shape.radii(theta, phi);
        if nargout == 1
          varargout{1} = [ r theta phi ];
        else
          varargout{1} = r;
          varargout{2} = theta;
          varargout{3} = phi;
        end
      end
    end
  end
end
