classdef EccentricSpheresNn < ott.drag.Stokes ...
    & ott.drag.mixin.VarViscosity ...
    & ott.drag.mixin.CalcInvDrag
% Calculate drag on an eccentric sphere using Gibson's NN approach.
% Inherits from :class:`Stokes`.
%
% Uses the NN from
%
%   Lachlan J. Gibson, et al. Phys. Rev. E 99, 043304
%   https://doi.org/10.1103/PhysRevE.99.043304
%
% Properties
%   - innerRadius   -- Radius of inner sphere
%   - outerRadius   -- Radius of outer sphere
%   - separation    -- Minimum separation between spheres
%   - viscosity     -- Viscosity of surrounding fluid (default: 1.0)
%   - forward       -- Computed drag tensor
%   - inverse       -- Computed from `forward`.
%
% See :class:`Stokes` for other methods/parameters.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    % Initial values are for construction/loading only (replaced on construct)
    innerRadius = 0     % Radius of inner sphere
    outerRadius = Inf   % Radius of outer sphere
    separation = 0      % Minimum separation between spheres
  end

  properties (Dependent)
    forwardInternal     % Drag on the inner sphere
  end

  methods
    function drag = EccentricSpheresNn(varargin)
      % Calculate the drag on an eccentric sphere.
      %
      % Usage:
      %   drag = EccentricSpheresNn(innerRadius, outerRadius, ...
      %     separation, viscosity, ...)
      %
      % Parameters:
      %   - innerRadius -- Radius of inner sphere
      %   - outerRadius -- Radius of outer sphere
      %   - separation -- Minimum separation between inner and outer sphere
      %   - viscosity -- Viscosity of medium (default: 1.0)
      %
      % Additional parameters are passed to corresponding class constructor.

      p = inputParser;
      p.addRequired('innerRadius', @isnumeric);
      p.addRequired('outerRadius', @isnumeric);
      p.addRequired('separation', @isnumeric);
      p.addParameter('viscosity', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      drag = drag@ott.drag.Stokes(unmatched{:});
      drag.innerRadius = p.Results.innerRadius;
      drag.outerRadius = p.Results.outerRadius;
      drag.separation = p.Results.separation;
      drag.viscosity = p.Results.viscosity;
    end
  end

  methods % Getters/setters
    function D = get.forwardInternal(drag)

      % Calculate stokes sphere drag
      D = ott.drag.StokesSphere(drag.innerRadius, drag.viscosity).forward;

      % Calculate parameters for Lachlan's script
      pD = drag.separation ./ drag.innerRadius;
      lam = drag.innerRadius ./ drag.outerRadius;

      [gy, fx, fxc, fz, gz] = WEES(pD, lam);

      % Add corrections for diagonal terms
      D(1, 1) = D(1, 1) .* fx;
      D(2, 2) = D(2, 2) .* fx;
      D(3, 3) = D(3, 3) .* fz;
      D(4, 4) = D(4, 4) .* gy;
      D(5, 5) = D(5, 5) .* gy;
      D(6, 6) = D(6, 6) .* gz;

      % Add cross-terms
      D(1, 5) = 6*pi*drag.viscosity*drag.innerRadius^2*fxc;
      D(2, 4) = -6*pi*drag.viscosity*drag.innerRadius^2*fxc;
      D(5, 1) = 8*pi*drag.viscosity*drag.innerRadius^2*fxc;
      D(4, 2) = -8*pi*drag.viscosity*drag.innerRadius^2*fxc;
    end

    function drag = set.innerRadius(drag, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'innerRadius must be positive numeric scalar');
      assert(val < drag.outerRadius, ...
          'innerRadius must be less than outerRadius');
      assert(drag.separation <= drag.outerRadius - val, ...
          'separation must be less that outerRadius-innerRadius');
      drag.innerRadius = val;
    end

    function drag = set.outerRadius(drag, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'outerRadius must be positive numeric scalar');
      assert(drag.innerRadius < val, ...
          'innerRadius must be less than outerRadius');
      assert(drag.separation <= val - drag.innerRadius, ...
          'separation must be less that outerRadius-innerRadius');
      drag.outerRadius = val;
    end

    function drag = set.separation(drag, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'separation must be positive numeric scalar');
      assert(val <= drag.outerRadius - drag.innerRadius, ...
          'separation must be less that outerRadius-innerRadius');
      drag.separation = val;
    end
  end
end

