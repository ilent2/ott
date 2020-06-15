classdef (Abstract) Extrusion < ott.shapes.Shape ...
    & ott.shapes.utils.CoordsCart
% Abstract class for shapes described by a height function.
% Inherits from :class:`utils.CoordsCart`, is aware of :class:`AxisymShape`.
%
% This class provides methods to help creating functions of the shapes
% where the height can be described as either::
%
%   z \equiv z(x, y)
%
% and::
%
%   z \equiv z(r)
%
% where :math:`r = \sqrt{x^2 + y^2}`.  If the shape inherits from
% :class:`AxisymShape` the second method is used.
%
% Abstract methods
%   - extrudeProfile   -- The extrusion profile
%   - get_maxXyRadius  -- Get the maximum radius in the XY plane
%
% Methods
%   - surf          -- Generate a visualisation of the shape
%   - normals       -- Calculate normals numerically
%   - insideXyz     -- Determines if a point is inside, calls extrudeProfile
%   - insideRtp     -- Spherical coordinate inputs, calls insideXyz
%   - get_maxRadius -- Numerically estimates maxRadius
%   - get_volume    -- Numerically estimates volume
%
% Properties
%   - maxRadius   -- maximum distance from shape origin
%   - maxXyRadius -- maximum radius in the xy-plane
%   - volume      -- volume of shape
%   - position    -- Location of shape ``[x, y, z]``

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    maxXyRadius  % maximum radius in the xy-plane
  end

  methods (Abstract)
    extrudeProfile    % The extrusion profile
    get_maxXyRadius   % Get the maximum radius in the XY plane
  end

  methods
    function varargout = surf(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Usage
      %   shape.surf(...) displays a visualisation of the shape
      %   in the current axes.
      %
      %   [X, Y, Z] = shape.surf() calculates the coordinates and arranges
      %   them in a grid suitable for use with matlab surf function.
      %
      % Optional named arguments
      %   - position (3x1 numeric) -- Offset for location of surface
      %     Default: ``[]``.
      %
      %   - rotation (3x3 numeric) -- Rotation matrix to apply to surface
      %     Default: ``[]``.
      %
      %   - npoints (2 numeric) -- Number of points to use for surface.
      %     If shape is a :class:`AxisymShape`, use ``[nr, nphi]``,
      %     else use ``[nx, ny]``.  Default: ``[50, 50]``.
      %
      %   - show_normals (logical) -- If true, plot normals at the surface
      %     locations using quiver.  Default: ``false``.
      %
      %   - axes (axes handle) -- axis to place surface in (default: gca)
      %
      %   - surfoptions {varargin} -- options to be passed to surf.

      p = inputParser;
      p.addParameter('npoints', [50, 50]);
      p.addParameter('surfoptions', {});
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('axes', []);
      p.addParameter('show_normals', false);
      p.parse(varargin{:});

      % Get the size from the user inputs
      sz = p.Results.npoints;
      if numel(sz) == 1
        sz = [sz sz];
      end

      % Generate coordinates
      if isa(shape, 'ott.shapes.AxisymShape')
        xrange = linspace(0, shape.maxXyRadius, sz(1));
        yrange = linspace(0, 2*pi, sz(2));
        [rr, rho] = meshgrid(xrange, yrange);

        xx = rr .* sin(rho);
        yy = rr .* cos(rho);
        zz = shape.extrudeProfile(rr);

      else
        xrange = linspace(-1, 1, sz(1)).*shape.maxXyRadius;
        yrange = linspace(-1, 1, sz(2)).*shape.maxXyRadius;
        [xx, yy] = meshgrid(xrange, yrange);
        zz = shape.extrudeProfile([xx(:), yy(:)].');

        % Output size changed: fix it
        if iscell(zz)
          zz{1} = reshape(zz{1}, size(xx));
          zz{2} = reshape(zz{2}, size(xx));
        else
          zz = reshape(zz, size(xx));
        end
      end

      % Add top and bottom of shape together
      xx = [xx, fliplr(xx)];
      yy = [yy, fliplr(yy)];
      if iscell(zz)
        zz = [zz{1}, fliplr(zz{2})];
      else
        zz = [zz, -fliplr(zz)];
      end

      % Draw the figure and handle rotations/translations
      [varargout{1:nargout}] = shape.surfCommon(p, sz, xx, yy, zz);
    end

    function nxyz = normalsXyzInternal(shape, xyz, varargin)
      % Calculate normals at the specified surface locations
      %
      % This method estimates the normal numerically.  If the normal
      % is easily computed, consider implementing a method in the
      % sub-class.
      %
      % Usage
      %   nxyz = shape.normalsXyz(xyz, ...) calculates the normal at
      %   the specified location.

      if isa(shape, 'ott.shapes.AxisymShape')

        dr = 1.0e-3 * shape.maxXyRadius;

        % Calculate normalised radial vectors (Cylindrical coordinates)
        vec_rr = xyz(1:2, :);
        rr = vecnorm(vec_rr);
        vec_rr = vec_rr ./ rr;      % For Cyl. to Cart.

        % Evaluate second location
        rr_dr = rr - dr;
        zshape = shape.extrudeProfile(rr_dr);

        if iscell(zshape)
          oldzshape = zshape{2};
          zshape = zshape{1};
          zshape(xyz(3, :) < 0) = oldzshape(xyz(3, :) < 0);
        else
          % Don't use sign, since original could be zero
          zshape(xyz(3, :) < 0) = -zshape(xyz(3, :) < 0);
        end

        % Calculate normalised surface normal (Cylindrical coordinates)
        vec_norm = zshape - xyz(3, :);
        vec_norm(2, :) = dr;
        vec_norm = vec_norm ./ vecnorm(vec_norm);
        vec_norm(:, xyz(3, :) < 0) = -vec_norm(:, xyz(3, :) < 0);

        % Calculate surface normal (Cartesian coordinates)
        nxyz = zeros(size(xyz));
        nxyz(3, :) = vec_norm(2, :);
        nxyz(1:2, :) = vec_norm(1, :) .* vec_rr;

      else

        dx = 1.0e-3 * shape.maxXyRadius;
        dy = 1.0e-3 * shape.maxXyRadius;

        xyz_dx = xyz + [dx; 0; 0];
        xyz_dy = xyz + [0; dy; 0];

        zshape_dx = shape.extrudeProfile(xyz_dx(1:2, :));
        zshape_dy = shape.extrudeProfile(xyz_dy(1:2, :));

        if iscell(zshape_dx)
          oldzshape = zshape_dx{2};
          zshape_dx = zshape_dx{1};
          zshape_dx(xyz(3, :) < 0) = oldzshape(xyz(3, :) < 0);

          oldzshape = zshape_dy{2};
          zshape_dy = zshape_dy{1};
          zshape_dy(xyz(3, :) < 0) = oldzshape(xyz(3, :) < 0);
        else
          % Don't use sign, since original could be zero
          zshape_dx(xyz(3, :) < 0) = -zshape_dx(xyz(3, :) < 0);
          zshape_dy(xyz(3, :) < 0) = -zshape_dy(xyz(3, :) < 0);
        end

        d1 = zeros(size(xyz));
        d2 = zeros(size(xyz));
        d1(3, :) = zshape_dx - xyz(3, :);
        d1(1, :) = dx;
        d2(3, :) = zshape_dy - xyz(3, :);
        d2(2, :) = dy;

        nxyz = cross(d1, d2);

      end
    end
  end

  methods % Getters/Setters
    function val = get.maxXyRadius(shape)
      val = shape.get_maxXyRadius();
    end
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz, varargin)
      % Determine if a point is within the shape

      % Get the shape profile
      zshape = shape.getExtrudeProfile(xyz(1:2, :));

      % Handle mirror symmetry
      if numel(zshape) == 2
        b = xyz(3, :) < zshape{2} & xyz(3, :) > zshape{1};
      else
        b = abs(xyz(3, :)) < zshape{1};
      end
    end

    function zshape = getExtrudeProfile(shape, xy)
      % Calculate shape profile
      %
      % This function calls extrudeProfile with one or two arguments
      % depending on the class type.  Returns a cell array.

      if isa(shape, 'ott.shapes.AxisymShape')
        rho = vecnorm(xy);
        zshape = shape.extrudeProfile(rho);

      else
        zshape = shape.extrudeProfile(xy);

      end

      if ~iscell(zshape)
        zshape = {zshape};
      end
    end

    function r = get_maxRadius(shape)
      % Estimate the maximum radius numerically (uses maxXyRadius)

      % TODO: Optimisation for Axisym shapes

      x = linspace(-shape.maxXyRadius, shape.maxXyRadius, 50);
      dx = diff(x(1:2));
      [xx, yy] = meshgrid(x, x);

      % Get the shape profile
      zshape = shape.getExtrudeProfile([xx(:), yy(:)].');

      if numel(zshape) == 2
        zheight = max(abs(zshape{1}), abs(zshape{2}));
      else
        zheight = abs(zshape{1});
      end

      zheight = reshape(zheight, size(xx));
      rr = sqrt(xx.^2 + yy.^2 + zheight.^2);

      r = max(rr(:)) + dx;
    end

    function v = get_volume(shape)
      % Estimate the volume numerically (uses maxXyRadius)

      % TODO: Optimisation for Axisym shapes

      x = linspace(-shape.maxXyRadius, shape.maxXyRadius, 50);
      dx = diff(x(1:2));
      [xx, yy] = meshgrid(x, x);

      % Get the shape profile
      zshape = shape.getExtrudeProfile([xx(:), yy(:)].');

      if numel(zshape) == 2
        zshape = zshape{2} - zshape{1};
      else
        zshape = 2.*zshape{1};
      end

      zshape = reshape(zshape, size(xx));
      zshape(isnan(zshape)) = 0;
      v = trapz(trapz(zshape)) * (x(2)-x(1)).^2;
    end
  end
end

