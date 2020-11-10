classdef AxisymFunc < ott.shape.Shape ...
    & ott.shape.mixin.AxisymShape ...
    & ott.shape.mixin.NumericalVolume ...
    & ott.shape.mixin.IntersectRayMarch
% Rotationally symmetry shape described by a function
%
% Properties
%   - func      -- Function describing surface
%   - type      -- Type of function (radial | angular | axial | axialSym)
%   - range     -- Range of function parameter values (default: [-Inf, Inf])
%
% Methods
%   - surf      -- Visualise the shape (via PatchMesh)
%
% Static methods
%   - BiconcaveDisc   -- Create a biconcave disk shape
%   - Pill            -- Create a pill tipped shaped cylinder
%
% Supported casts
%   - AxisymInterp
%   - PatchMesh       -- Via AxisymInterp

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    func        % Function describing surface
    type        % Type of function (radial | angular | axial | axialSym)
  end
  
  properties (Hidden, SetAccess=protected)
    rangeInternal     % Internal property set by range
  end

  properties (Dependent)
    range       % Range of function parameter values (default: [-Inf, Inf])
    boundingBox       % Box surrounding shape
    perimeter         % Perimeter of shape
    maxRadius         % Maximum radius of the particle
    starShaped        % True if the particle is star-shaped
    xySymmetry        % True if the particle is xy-plane mirror symmetric
  end

  methods (Static)
    function shape = BiconcaveDisc(varargin)
      % Construct a biconcave disc shape
      %
      % This shape can be used to model cells such as unstressed
      % Red Blood Cells.  It implements the function::
      %
      %   z(r) = D \sqrt{1 - \frac{4r^2}{D^2}} \left(a_0 +
      %       \frac{a_1 r^2}{D^2} + \frac{a_2 r^4}{D^4} \right)
      %
      % where :math:`D` is the particle diameter and :math:`a` are shape
      % coefficients.
      %
      % Usage
      %   shape = AxisymFunc.BiconcaveDisc(radius, coefficients, ...)
      %
      % Parameters
      %   - radius (numeric) -- Radius of disc.  Default: ``7.82``.
      %   - coefficients (3 numeric) -- Coefficients describing shape.
      %     Default: ``[0.0518, 2.0026, -4.491]``.
      %
      % Additional parameters are passed to shape constructor.

      p = inputParser;
      p.addOptional('radius', 7.82);
      p.addOptional('coefficients', [0.0518, 2.0026, -4.491]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      assert(numel(p.Results.coefficients) == 3, ...
          'coefficients must be 3 element numeric vector');
      a0 = p.Results.coefficients(1);
      a1 = p.Results.coefficients(2);
      a2 = p.Results.coefficients(3);

      assert(isscalar(p.Results.radius), ...
          'radius must be numeric scalar');
      radius = p.Results.radius;

      func = @(r) 2.0 .* radius .* sqrt(1 - r.^2./radius.^2) ...
        .* (a0 + (a1.*r.^2) ./ radius.^2 ./ 4 ...
        + (a2.*r.^4) ./ radius.^4 ./ 16);

      shape = ott.shape.AxisymFunc('func', func, 'type', 'axialSym', ...
          'range', [0, radius], unmatched{:});
    end

    function shape = Pill(varargin)
      % Construct a pill-shaped particle
      %
      % Constructs a cylindrical shaped rod with spherical end caps.
      %
      % Usage
      %   shape = AxisymFunc.Pill(height, radius, capRadius, ...)
      %
      % Parameters
      %   - height (numeric) -- Total height of pill including
      %     end caps.  Default: ``2.0``.
      %
      %   - radius (numeric) -- Radius of rod segment.  Default: ``0.5``.
      %
      %   - capRadius (numeric) -- Radius of end cap.  Default: ``radius``.
      %     If `capRadius` is greater than `radius`, adds a sharp edge
      %     at the intersection of the cap and the rod.  If radius is
      %     less, adds a flat end-cap.
      %
      % Additional parameters are passed to class constructor.

      p = inputParser;
      p.addOptional('height', 2.0);
      p.addOptional('radius', 0.5);
      p.addOptional('capRadius', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      height = p.Results.height;
      radius = p.Results.radius;
      capRadius = p.Results.capRadius;
      if isempty(capRadius)
        capRadius = radius;
      end

      assert(isscalar(height) && isnumeric(height) && height > 0, ...
          'height must be numeric scalar > 0');
      assert(isscalar(radius) && isnumeric(radius) && radius > 0, ...
          'radius must be numeric scalar > 0');
      assert(isscalar(capRadius) && isnumeric(capRadius) && capRadius > 0, ...
          'capRadius must be numeric scalar > 0');

      func = @(r) height/2 - capRadius ...
          + sqrt(capRadius.^2 - (max(0, r - max(0, radius - capRadius))).^2);
      shape = ott.shape.AxisymFunc('func', func, 'type', 'axialSym', ...
          'range', [0, radius], unmatched{:});
    end
  end

  methods
    function shape = AxisymFunc(varargin)
      % Construct a new rotationally symmetric shape from a function
      %
      % Usage
      %   shape = AxisymFunc(func, type, ...)
      %
      % Parameters
      %   - func (function_handle) -- A function handle describing
      %     the shape surface.  The function should take a single
      %     vectorised argument.  The argument will depend on the
      %     `type`: radial: z, angular: theta, axial: r
      %
      %   - type (enum) -- Type of function.  Can either be
      %     'angular', 'radial', 'axial' or 'axialSym'.
      %
      % Optional named arguments
      %   - range (2 numeric) -- Range of function parameter values.
      %     Default: ``[-Inf, Inf]`` (radial), ``[-pi, pi]`` (angular),
      %     and ``[0, Inf]`` (axial/axialSym).
      %
      % Additional parameters passed to base.

      p = inputParser;
      p.addOptional('func', [], @(x) isa(x, 'function_handle'));
      p.addOptional('type', [], @ischar);
      p.addParameter('range', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = shape@ott.shape.Shape(unmatched{:});
      shape.func = p.Results.func;
      shape.type = p.Results.type;

      if isempty(p.Results.range)
        if strcmpi(shape.type, 'angular')
          shape.range = [-pi, pi];
        elseif strcmpi(shape.type, 'radial')
          shape.range = [-Inf, Inf];
        else
          shape.range = [0, Inf];
        end
      else
        shape.range = p.Results.range;
      end
    end

    function varargout = surfPoints(shape, varargin)
      % Cast to PatchMesh and call surfPoints
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)
      %
      % Optional named parameters
      %   - resolution ([nfunc, nphi]) -- Number of faces in the function
      %     direction and around the axis.  Default: ``[20, 20]``.
      %
      % Additional named parameters are passed to PatchMesh.surf.

      p = inputParser;
      p.addParameter('resolution', [20, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shape.PatchMesh(shape, 'resolution', p.Results.resolution);
      [varargout{1:nargout}] = shape.surfPoints(unmatched{:});
    end

    function shape = ott.shape.AxisymInterp(shape, varargin)
      % Cast shape to a AxisymInterp
      %
      % Usage
      %   shape = ott.shape.AxisymInterp(shape, ...)
      %
      % Optional named parameters
      %   - resolution (numeric) -- Number of faces in the function
      %     direction. Default: ``20``.
      %
      % Additional named parameters are passed to constructor.

      p = inputParser;
      p.addParameter('resolution', 20);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Evaluate function
      A = linspace(shape.range(1), shape.range(2), p.Results.resolution+1);
      V = shape.func(A);

      % Convert to cylindrical coordinates
      switch shape.type
        case 'angular'
          points = [V.*cos(A); V.*sin(A)];
        case 'radial'
          points = [V; A];
        case 'axial'
          points = [A; V];
        case 'axialSym'
          points = [A, fliplr(A); V, -fliplr(V)];

          % Remove duplicates
          points = unique(points.', 'rows', 'stable').';
        otherwise
          error('Unknown shape type');
      end

      % Construct AxisymInterp
      shape = ott.shape.AxisymInterp(points, ...
          'position', shape.position, 'rotation', shape.rotation, ...
          unmatched{:});
    end

    function shape = ott.shape.PatchMesh(shape, varargin)
      % Cast shape to a PatchMesh via AxisymInterp
      %
      % Usage
      %   shape = ott.shape.PatchMesh(shape, ...)
      %
      % Optional named parameters
      %   - resolution ([nfunc, nphi]) -- Number of faces in the function
      %     direction and around the axis.  Default: ``[20, 20]``.
      %
      % Additional named parameters are passed to the PatchMesh cast.

      p = inputParser;
      p.addParameter('resolution', [20, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Cast to AxisymInterp
      shape = ott.shape.AxisymInterp(shape, ...
          'resolution', p.Results.resolution(1));

      % Cast to patch
      shape = ott.shape.PatchMesh(shape, ...
          'resolution', p.Results.resolution(2), unmatched{:});
    end
  end

  methods (Hidden)
    function b = insideRtInternal(shape, rt)
      switch shape.type
        case 'angular'
          mask = rt(2, :) >= shape.range(1) & rt(2, :) <= shape.range(2);
          b = mask;
          b(mask) = rt(1, mask) <= shape.func(rt(2, mask));
        otherwise
          rz = rt(1, :).*[cos(rt(2, :)); sin(rt(2, :))];
          b = shape.insideRzInternal(rz);
      end
    end

    function b = insideRzInternal(shape, rz)
      switch shape.type
        case 'angular'
          rt = [vecnorm(rz); atan2(rz(1, :), rz(2, :))];
          b = shape.insideRtInternal(rt);
        case 'radial'
          mask = rz(2, :) >= shape.range(1) & rz(2, :) <= shape.range(2);
          b = mask;
          b(mask) = rz(1, mask) <= shape.func(rz(2, mask));
        case 'axial'
          mask = rz(1, :) >= shape.range(1) & rz(1, :) <= shape.range(2);
          b = mask;
          b(mask) = rz(2, mask) <= shape.func(rz(1, mask));
        case 'axialSym'
          mask = rz(1, :) >= shape.range(1) & rz(1, :) <= shape.range(2);
          b = mask;
          b(mask) = abs(rz(2, mask)) <= shape.func(rz(1, mask));
        otherwise
          error('Unknown shape type');
      end
    end

    function nxz = normalsRtInternal(shape, rt)
      switch shape.type
        case 'angular'
          mask = rt(2, :) >= shape.range(1) & rt(2, :) <= shape.range(2);

          dt = 1.0e-3;

          % Use Backward difference for points near end
          use_bk = rt(2, :) >= shape.range(2)-dt;
          rt(2, use_bk) = rt(2, use_bk) - dt;

          nxz = zeros(2, size(rt, 2));
          nxz(:, ~mask) = nan;

          theta = rt(2, mask);
          tp = (shape.func(theta+dt) - shape.func(theta));
          nxz(1, mask) = tp.*cos(theta) + dt.*sin(theta);
          nxz(2, mask) = tp.*sin(theta) - dt.*cos(theta);

          nxz(:, mask) = nxz(:, mask) ./ vecnorm(nxz(:, mask));
        otherwise
          rz = rt(1, :).*[cos(rt(2, :)); sin(rt(2, :))];
          nxz = shape.normalsRzInternal(rz);
      end
    end

    function nxz = normalsRzInternal(shape, rz)
      % Estimate normal numerically

      switch shape.type
        case 'angular'
          rt = [vecnorm(rz); atan2(rz(1, :), rz(2, :))];
          nxz = shape.normalsRtInternal(rt);
        case 'radial'
          mask = rz(2, :) >= shape.range(1) & rz(2, :) <= shape.range(2);

          dz = 1.0e-3;

          % Use Backward difference for points near end
          use_bk = rz(2, :) >= shape.range(2)-dz;
          rz(2, use_bk) = rz(2, use_bk) - dz;

          nxz = zeros(2, size(rz, 2));
          nxz(:, ~mask) = nan;
          nxz(1, mask) = dz;
          nxz(2, mask) = -(shape.func(rz(2, mask)+dz) ...
              - shape.func(rz(2, mask)));
          nxz(:, mask) = nxz(:, mask) ./ vecnorm(nxz(:, mask));

        case 'axial'
          mask = rz(1, :) >= shape.range(1) & rz(1, :) <= shape.range(2);

          dr = 1.0e-3;

          % Use Backward difference for points near end
          use_bk = rz(1, :) >= shape.range(2)-dr;
          rz(1, use_bk) = rz(1, use_bk) - dr;

          nxz = zeros(2, size(rz, 2));
          nxz(:, ~mask) = nan;
          nxz(1, mask) = -(shape.func(rz(1, mask)+dr) ...
              - shape.func(rz(1, mask)));
          nxz(2, mask) = dr;
          nxz(:, mask) = nxz(:, mask) ./ vecnorm(nxz(:, mask));

          nxz(:, rz(2, :) <= 0) = [0; -1];

        case 'axialSym'
          mask = rz(1, :) >= shape.range(1) & rz(1, :) <= shape.range(2);

          dr = 1.0e-3;

          % Use Backward difference for points near end
          use_bk = rz(1, :) >= shape.range(2)-dr;
          rz(1, use_bk) = rz(1, use_bk) - dr;

          nxz = zeros(2, size(rz, 2));
          nxz(:, ~mask) = nan;
          nxz(1, mask) = -(shape.func(rz(1, mask)+dr) ...
              - shape.func(rz(1, mask)));
          nxz(2, mask) = dr;
          nxz(:, mask) = nxz(:, mask) ./ vecnorm(nxz(:, mask));

          nxz(2, rz(2, :) <= 0) = -nxz(2, rz(2, :) <= 0);

        otherwise
          error('Unknown shape type');
      end
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
      %   - resolution ([nfunc, nphi]) -- Number of faces in the function
      %     direction and around the axis.  Default: ``[20, 20]``.
      %
      % Additional named parameters are passed to PatchMesh.surf.

      p = inputParser;
      p.addParameter('resolution', [20, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shape.PatchMesh(shape, 'resolution', p.Results.resolution);
      S = shape.surfInternal(unmatched{:});
    end

    function shape = scaleInternal(shape, sc)
      if strcmpi(shape.type, 'angular')
        shape.func = @(x) shape.func(x) .* sc;
      else
        shape.range = shape.range * sc;
        shape.func = @(x) shape.func(x ./ sc) .* sc;
      end
    end
  end

  methods % Getters/setters
    function shape = set.func(shape, val)
      assert(isa(val, 'function_handle'), ...
          'func should be a function_handle');
      shape.func = val;
    end

    function shape = set.type(shape, val)
      assert(any(strcmpi(val, {'radial', 'axial', 'axialSym', 'angular'})), ...
          'type should be ''radial''|''axial''|''axialSym''|''angular''');
      shape.type = val;
    end

    function shape = set.range(shape, val)
      assert(isnumeric(val) && isvector(val) && numel(val) == 2, ...
          'range must be 2 element numeric vector');
      assert(val(1) < val(2), 'range must be sorted in ascending order');

      % Type specific checks
      if strcmpi(shape.type, 'axial') || strcmpi(shape.type, 'axialSym')
        assert(val(1) >= 0, 'range(1) must be >= 0 for axial types');
      elseif strcmpi(shape.type, 'angular')
        assert(val(1) >= -pi && val(2) <= pi, ...
            'range must be in [-pi, pi] for angular type');
      end

      shape.rangeInternal = val(:).';
    end
    function val = get.range(shape)
      val = shape.rangeInternal;
    end

    function bb = get.boundingBox(shape)
      switch shape.type
        case 'radial'
          [~, R] = fminbnd(@(x) -shape.func(x), ...
              shape.range(1), shape.range(2));
          bb = [R, -R; R, -R; shape.range(1:2)];
        case 'axial'
          [~, Z] = fminbnd(@(x) -shape.func(x), ...
              shape.range(1), shape.range(2));
          R = shape.range(2);
          bb = [-R, R; -R, R; 0, -Z];
        case 'axialSym'
          [~, Z] = fminbnd(@(x) -shape.func(x), ...
              shape.range(1), shape.range(2));
          R = shape.range(2);
          bb = [-R, R; -R, R; Z, -Z];
        case 'angular'
          [~, R] = fminbnd(@(x) -shape.func(x).*cos(x), ...
              shape.range(1), shape.range(2));
          [~, pZ] = fminbnd(@(x) -shape.func(x).*sin(x), ...
              shape.range(1), shape.range(2));
          [~, mZ] = fminbnd(@(x) shape.func(x).*sin(x), ...
              shape.range(1), shape.range(2));
          bb = -[-R, R; -R, R; mZ, pZ];
        otherwise
          error('Unknown shape type');
      end
    end

    function r = get.maxRadius(shape)
      if strcmpi(shape.type, 'angular')
        [~, r] = fminbnd(@(x) -shape.func(x), ...
            shape.range(1), shape.range(2));
      else
        [~, r] = fminbnd(@(x) -sqrt(shape.func(x).^2 + x.^2), ...
            shape.range(1), shape.range(2));
      end
    end

    function b = get.starShaped(shape)
      switch shape.type
        case 'angular'
          b = true;
        case 'radial'
          % Check the derivative
          dt = 1.0e-3;
          [~, grad] = fminbnd(@(x) atan2(x+dt, shape.func(x+dt)) ...
              - atan2(x, shape.func(x)), ...
              shape.range(1), shape.range(2)-dt);
          b = grad > 0;
        case 'axial'
          % Check derivative, assumes f(x) > 0 V x
          dt = 1.0e-3;
          [~, grad] = fminbnd(@(x) atan2(shape.func(x+dt), x+dt) ...
              - atan2(shape.func(x), x), ...
              shape.range(1), shape.range(2)-dt);
          b = grad > 0;
        case 'axialSym'
          % Check derivative, assumes f(x) > 0 V x
          dt = 1.0e-3;
          [~, grad] = fminbnd(@(x) atan2(shape.func(x+dt), x+dt) ...
              - atan2(shape.func(x), x), ...
              shape.range(1), shape.range(2)-dt);
          b = grad > 0;
        otherwise
          error('Unknown shape type');
      end
    end

    function b = get.xySymmetry(shape)
      switch shape.type
        case 'angular'
          b = shape.range(1) == -shape.range(2);
          if b
            [~, v] = fminbnd(@(x) shape.func(x) - shape.func(-x), ...
                0, shape.range(2));
            tol = 1.0e-8;
            b = abs(v) <= tol;
          end
        case 'radial'
          b = shape.range(1) == -shape.range(2);
          if b
            [~, v] = fminbnd(@(x) shape.func(x) - shape.func(-x), ...
                0, shape.range(2));
            tol = 1.0e-8;
            b = abs(v) <= tol;
          end
        case 'axialSym'
          b = true;
        case 'axial'
          b = false;
        otherwise
          error('Unknown shape type');
      end
    end
  end
end

