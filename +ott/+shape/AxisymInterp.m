classdef AxisymInterp < ott.shape.Shape ...
    & ott.shape.mixin.AxisymShape ...
    & ott.shape.mixin.NumericalVolume ...
    & ott.shape.mixin.IntersectRayMarch
% Rotationally symmetric shape described by discrete set of points.
%
% Shape produced when converting this object to a patch uses linear
% interpolation between points and discrete rotational segments.
%
% Properties
%   - points        -- Discrete points describing surface [rho; z]
%
% Methods
%   - boundaryPoints  -- Calculate points for line integration
%   - surfPoints      -- Calculate points for surface integration
%
% Supported casts
%   - PatchMesh
%
% Static methods
%   - Bicone             -- Create a bicone
%   - ConeTippedCylinder -- Create a cone-tipped cylinder
%
% Volume is computed numerically, may change in future.
%
% Additional properties/methods inherited from base.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    points          % Discrete points describing surface [rho; z]
  end

  properties (Dependent)
    maxRadius          % Maximum particle radius
    perimeter          % Perimeter of shape
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    starShaped         % True if the particle is star-shaped
    xySymmetry         % True if the particle is xy-plane mirror symmetric
  end

  methods (Static)
    function shape = ConeTippedCylinder(varargin)
      % Construct a cone-tipped cylinder
      %
      % Usage
      %   shape = ConeTippedCylinder(height, radius, coneHeight, ...)
      %
      % Parameters
      %   - height (numeric) -- total height of the shape (default: 2.0)
      %   - radius (numeric) -- radius of cylinder (default: 1.0)
      %   - coneHeight (numeric) -- Height of the cone segment (default: 0.5)
      %
      % Additional parameters passed to class constructor.

      p = inputParser;
      p.addOptional('height', 2.0, @isnumeric);
      p.addOptional('radius', 1.0, @isnumeric);
      p.addOptional('coneHeight', 0.5, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      radius = p.Results.radius;
      coneHeight = p.Results.coneHeight;
      height = p.Results.height;

      assert(isscalar(radius) && isnumeric(radius) && radius > 0, ...
          'radius must be positive numeric scalar');
      assert(isscalar(height) && isnumeric(height) && height > 0, ...
          'height must be positive numeric scalar');
      assert(isscalar(coneHeight) && isnumeric(coneHeight) ...
          && coneHeight > 0 && coneHeight < height/2, ...
          'coneHeight must be numeric scalar in range (0, height/2)');

      points = [0, radius, radius, 0; ...
          height/2, height/2-coneHeight, -height/2+coneHeight, -height/2];
      shape = ott.shape.AxisymInterp(points, unmatched{:});
    end

    function shape = Bicone(varargin)
      % Construct a bicone
      %
      % A Bicone is xy-mirror symmetric and has three vertices,
      % two on the +(ve)/-(ve) axes and one in the mirror symmetric plane.
      %
      % Usage
      %   shape = AxisymInterp.Bicone(height, radius, ...)
      %
      % Parameters
      %   - height (numeric) -- Total height of the shape (default: 2.0)
      %   - radius (numeric) -- Radius of the cone (default: 1.0)
      %
      % Additional parameters passed to class constructor.

      p = inputParser;
      p.addOptional('height', 2.0, @isnumeric);
      p.addOptional('radius', 1.0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      radius = p.Results.radius;
      height = p.Results.height;

      assert(isscalar(radius) && isnumeric(radius) && radius > 0, ...
          'radius must be positive numeric scalar');
      assert(isscalar(height) && isnumeric(height) && height > 0, ...
          'height must be positive numeric scalar');

      points = [0, radius, 0; height/2, 0, -height/2];
      shape = ott.shape.AxisymInterp(points, unmatched{:});
    end
  end

  methods
    function shape = AxisymInterp(varargin)
      % Construct a new rotationally symmetry shape from discrete points
      %
      % Usage
      %   shape = AxisymInterp(points, ...)
      %
      % Parameters
      %   - points (2xN numeric) -- Array of points in cylindrical coordinates
      %     describing shape geometry. [rho; z]
      %
      % Additional parameters passed to base.

      p = inputParser;
      p.addOptional('points', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shape.Shape(unmatched{:});
      shape.points = p.Results.points;
    end

    function [xyz, nxyz, ds] = boundaryPoints(shape, npts)
      % Calculate points for line integration
      %
      % Usage
      %   [xyz, nxyz, ds] = shape.boundaryPoints(npts, ...)
      %
      % Parameters
      %   - npts (numeric) -- Number of equally spaced points along surface.

      % Calculate the point spacing
      ds = shape.perimeter / npts / 2.0;

      % The following is based on axisym_boundarypoints from OTTv1

      % Calculate length of each line segment
      s=sqrt((rho(2:end)-rho(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2);

      %Don't ask me how this works. It does. It's simple algebra in the end...
      %and yes, it can be done more intelligently.
      zout=zeros(npts,1);
      rhoout=zout;
      nxyz=zeros(3, npts);

      sdeficit=0;
      ncum=0;

      for ii=2:length(rho)
          N=s(ii-1)/ds;
          Nused=round(N+sdeficit);

          nc=[-(z(ii)-z(ii-1));0;(rho(ii)-rho(ii-1))];
          nc=nc./vecnorm(nc);

          if Nused>=1
              drho=(rho(ii)-rho(ii-1))/N*ones(Nused,1);

              rhot=cumsum(drho)-drho/2-sdeficit*drho(1);
              rhoout(ncum+(1:Nused))=rho(ii-1)+rhot;

              dz=(z(ii)-z(ii-1))/N*ones(Nused,1);

              zt=cumsum(dz)-dz/2-sdeficit*dz(1);
              zout(ncum+(1:Nused))=z(ii-1)+zt;

              nxyz(:, ncum+(1:Nused))=repmat(nc,[length(zt),1]);

              sdeficit=(N-Nused+sdeficit);
          else
              sdeficit=sdeficit+N;
          end

          ncum=ncum+Nused;

      end

      % Truncate points if not allocated
      if ncum < npts
        warning('ott:shapes:AxisymInterp:boundarypoints_length', ...
            'Number of points generated does not match request');
        zout = zout(1:ncum);
        rhoout = zout(1:ncum);
        nxyz = nxyz(:, 1:ncum);
      end

      % Package output
      xyz = [rhoout,zeros(size(rhoout)),zout].';

      % Calculate area
      if nargout >= 3
        dst=zeros(size(rtp, 1),3);

        %calcultes area elements
        dst(2:end-1,1)=(rhoout(3:end)-rhoout(1:end-2))/2;
        dst(2:end-1,3)=(zout(3:end)-zout(1:end-2))/2;
        dst(1,1)=(mean(rhoout(1:2))-rho(1));
        dst(1,3)=(mean(zout(1:2))-z(1));
        dst(end,1)=(rho(end)-mean(rhoout(end-1:end)));
        dst(end,3)=(z(end)-mean(zout(end-1:end)));

        % a general axisymmetric conic region has the
        % following area (sans factor of 2*pi):
        ds=(rtp(:,1).*sqrt(sum(abs(dst).^2,2)).*sin(rtp(:,2)));
      end
    end

    function varargout = surfPoints(shape, varargin)
      % Cast to PatchMesh and call surfPoints
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)
      %
      % Optional named parameters
      %   - resolution (numeric) -- Number of faces in angular direction.
      %     Default: ``20``.
      %
      % Additional named parameters are passed to PatchMesh.surf.

      p = inputParser;
      p.addParameter('resolution', 20);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shape.PatchMesh(shape, 'resolution', p.Results.resolution);
      [varargout{1:nargout}] = shape.surfPoints(unmatched{:});
    end

    function shape = ott.shape.PatchMesh(shape, varargin)
      % Cast shape to a PatchMesh
      %
      % Usage
      %   shape = ott.shape.PatchMesh(shape, ...)
      %
      % Optional named parameters
      %   - resolution (numeric) -- Number of faces in rotationally
      %     symmetric directions.  Default: ``20``.
      %
      % Additional named parameters are passed to constructor.

      p = inputParser;
      p.addParameter('resolution', 20);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      phi = linspace(0, 2*pi, p.Results.resolution+1).';
      X = shape.points(1, :).*sin(phi);
      Y = shape.points(1, :).*cos(phi);
      Z = repmat(shape.points(2, :), size(phi, 1), 1);

      shape = ott.shape.PatchMesh.FromSurfMatrix(X, Y, Z, ...
          'position', shape.position, 'rotation', shape.rotation, ...
          unmatched{:});
    end
  end

  methods (Hidden)
    function b = insideRtInternal(shape, rt)
      % Use insideRzInternal instead
      rz = rt(1, :).*[sin(rt(2, :)), cos(rt(2, :))];
      b = shape.insideRzInternal(rz);
    end

    function b = insideRzInternal(shape, rz)
      % Determine if inside by counting how many edges are above point
      % "up" is in the radial direction.

      % Find points greater or less than query points
      %   - ignore points with same z position
      %   - only care about points connected to edges
      pts_lt = shape.points(2, 1:end-1) < rz(2, :).';
      pts_gt = shape.points(2, 2:end) > rz(2, :).';

      % Find edges crossing points
      edges = pts_lt & pts_gt;

      % Calculate intersection for every edge
      m = (shape.points(1, 2:end) - shape.points(1, 1:end-1)) ...
          ./ (shape.points(2, 2:end) - shape.points(2, 1:end-1));
      x = rz(2, :).' - shape.points(2, 1:end-1);
      c = shape.points(1, 1:end-1);
      y4 = m.*x + c;

      % Discard edges we don't care about
      y4 = y4 .* edges;
      y4(y4 < rz(1, :)) = 0;

      % Determine if inside (count edges above)
      b = mod(sum(y4 > 0, 2).', 2) == 1;
    end

    function nz = normalsRtInternal(shape, rt)
      % Defer to normalsRzInternal
      rz = rt(1, :).*[sin(rt(2, :)), cos(rt(2, :))];
      nz = shape.normalsRzInternal(rz);
    end

    function nz = normalsRzInternal(shape, rz)
      % Calculate normals by finding the closest edge to the point
      % This follows a similar procedure to the insideRzInternal function.

      % Find points greater or less than query points
      %   - ignore points with same z position
      %   - only care about points connected to edges
      pts_lt = shape.points(2, 1:end-1) < rz(2, :).';
      pts_gt = shape.points(2, 2:end) > rz(2, :).';

      % Find edges crossing points
      edges = pts_lt & pts_gt;

      % Calculate intersection for every edge
      m = (shape.points(1, 2:end) - shape.points(1, 1:end-1)) ...
          ./ (shape.points(2, 2:end) - shape.points(2, 1:end-1));
      x = rz(2, :).' - shape.points(2, 1:end-1);
      c = shape.points(1, 1:end-1);
      y4 = m.*x + c;

      % Discard edges we don't care about
      y4 = y4 .* edges;

      % Find which edge is closest
      [~, idx] = min(abs(y4 - rz(1, :)), [], 2);

      % Compute formals and normalize
      nz = [m(idx).'; ones(size(idx)).'];
      nz = nz ./ vecnorm(nz);
    end

    function S = surfInternal(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Converts the shape to a PatchMesh and calls surf.
      %
      % Usage
      %   S = shape.surfInternal(...)
      %
      % Optional named parameters
      %   - resolution (numeric) -- Number of faces in angular direction.
      %     Default: ``20``.
      %
      % Additional named parameters are passed to PatchMesh.surfInternal.

      p = inputParser;
      p.addParameter('resolution', 20);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shape.PatchMesh(shape, 'resolution', p.Results.resolution);
      S = shape.surfInternal(unmatched{:});
    end

    function shape = scaleInternal(shape, sc)
      shape.points = shape.points .* sc;
    end
  end

  methods % Getters/setters
    function shape = set.points(shape, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 2, ...
        'points must be numeric 2xN matrix');
      assert(all(val(1, :) >= 0), 'points(1, :) must all be positive');
      shape.points = val;
    end

    function r = get.maxRadius(shape)
      r = max(vecnorm(shape.points));
    end

    function p = get.perimeter(shape)
      p = sum(vecnorm(diff(shape.points, [], 2)));
    end

    function bb = get.boundingBox(shape)
      X = max(shape.points(1, :));
      bb = [-X, X; -X, X; min(shape.points(2, :)), max(shape.points(2, :))];
    end

    function b = get.starShaped(shape)
      angles = atan2(shape.points(2, :), shape.points(1, :));
      b = issorted(angles, 'monotonic');
    end

    function b = get.xySymmetry(shape)
      pts = shape.points;
      npts = fliplr([shape.points(1, :); -shape.points(2, :)]);
      b = all(all(pts == npts));
    end
  end
end
