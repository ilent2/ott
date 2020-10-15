classdef StokesSphere < ott.drag.Stokes ...
    & ott.drag.mixin.CalcInvDrag ...
    & ott.drag.mixin.VarViscosity
% Drag tensor for a sphere with Stokes Drag
%
% Properties
%   - radius    -- Radius of sphere
%   - viscosity -- Viscosity of medium
%   - forward   -- Calculated drag tensor
%   - inverse   -- Calculate from forward
%
% Supported casts
%   - ott.shape.Shape -- Constructs a sphere shape
%
% See :class:`Stokes` for other methods/parameters.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius        % Radius of sphere
  end

  properties (Dependent)
    forwardInternal       % Calculated drag tensor
  end

  methods
    function drag = StokesSphere(varargin)
      % Calculate drag tensors for spherical particle in Stokes drag.
      %
      % Usage:
      %   tensor = Sphere(radius, viscosity, ...)
      %
      % Parameters
      %   - radius    -- (numeric) Radius of particle (default: 1.0)
      %   - viscosity -- (numeric) Viscosity of medium (default: 1.0)
      %
      % Parameters can also be passed as named arguments.
      % Unmatched parameters are passed to :class:`Stokes`.

      p = inputParser;
      p.addOptional('radius', 1.0, @isnumeric);
      p.addOptional('viscosity', 1.0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      drag = drag@ott.drag.Stokes(unmatched{:});
      drag.viscosity = p.Results.viscosity;
      drag.radius = p.Results.radius;
    end

    function shape = ott.shape.Shape(drag)
      % Construct a geometrical representation of the shape
      shape = ott.shape.Sphere(drag.radius, 'rotation', drag.rotation);
    end
  end

  methods % Getters/setters
    function D = get.forwardInternal(drag)
      D = diag([6*pi*drag.viscosity*drag.radius.*[1;1;1]; ...
                8*pi*drag.viscosity*drag.radius.^3.*[1;1;1]]);
    end

    function drag = set.radius(drag, val)
      assert(isnumeric(val) && isscalar(val), ...
          'radius must be numeric scalar');
      drag.radius = val;
    end
  end
end
