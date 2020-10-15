classdef StokesCylinder < ott.drag.Stokes ...
    & ott.drag.mixin.CalcInvDrag ...
    & ott.drag.mixin.VarViscosity
% Drag tensor for long/slender cylindrical particles.
%
% Uses the results from
%
%   Maria M. Tirado and José García de la Torre
%   J. Chem. Phys. 71, 2581 (1979); https://doi.org/10.1063/1.438613
%
% And from
%
%   María M. Tirado and José García de la Torre
%   J. Chem. Phys. 73, 1986 (1980); https://doi.org/10.1063/1.440288
%
% which provide lookup tables for the drag corrections for slender
% cylinders (aspect ratio `height/diameter` above 2).
%
% Properties
%   - radius      -- Cylinder radius
%   - height      -- Cylinder height
%   - viscosity   -- Medium viscosity
%   - forward     -- Calculated drag tensor
%   - inverse     -- Inverse drag tensor (calculated from `forward`)
%
% Supported casts
%   - ott.shape.Shape -- Construct a cylinder
%
% Additional properties/methods inherited from :class:`Stokes`.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius      % Cylinder radius (numeric)
    height      % Cylinder height (numeric)
  end

  properties (Dependent)
    forwardInternal     % Forward drag tensor
  end

  properties (Constant, Hidden)
    
    % Axial rotation scaling factor (from paper)
    A0 = 3.841;

    % Rotational drag coefficients (table 2 from 10.1063/1.440288)
    rotationTable = [
        0.50, -0.216, 0.294;
        0.45, -0.260, 0.269;
        0.40, -0.303, 0.243;
        0.35, -0.348, 0.216;
        0.30, -0.392, 0.189;
        0.25, -0.436, 0.159;
        0.20, -0.481, 0.130;
        0.15, -0.526, 0.099;
        0.10, -0.571, 0.067;
        0.05, -0.616, 0.034;
        0, -0.662, 0];

    % Translational drag coefficients (table 1 from 10.1063/1.438613)
    translationTable = [
        0.50, 0.99, 0.25, 0.62;
        0.46, 0.97, 0.21, 0.59;
        0.42, 0.96, 0.18, 0.57;
        0.38, 0.94, 0.14, 0.54;
        0.33, 0.92, 0.10, 0.51;
        0.30, 0.91, 0.07, 0.49;
        0.25, 0.90, 0.03, 0.46;
        0.20, 0.88, -0.02, 0.43;
        0.17, 0.87, -0.05, 0.41;
        0.10, 0.86, -0.11, 0.37;
        0.06, 0.85, -0.15, 0.35;
        0.02, 0.84, -0.18, 0.33;
        0.00, 0.84, -0.20, 0.32];
  end

  methods
    function drag = StokesCylinder(varargin)
      % Calculate drag tensors for cylindrical particles
      %
      % Usage:
      %   tensor = Sphere(radius, height, viscosity, ...)
      %
      % Parameters
      %   - radius    -- (numeric) Radius of cylinder (default: 1.0)
      %   - height    -- (numeric) Height of cylinder (default: 4.0)
      %   - viscosity -- (numeric) Viscosity of medium (default: 1.0)
      %
      % Parameters can also be passed as named arguments.
      % Unmatched parameters are passed to :class:`Stokes`.

      p = inputParser;
      p.addOptional('radius', 1.0, @isnumeric);
      p.addOptional('height', 1.0, @isnumeric);
      p.addOptional('viscosity', 1.0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      drag = drag@ott.drag.Stokes(unmatched{:});
      drag.viscosity = p.Results.viscosity;
      drag.radius = p.Results.radius;
      drag.height = p.Results.height;
    end

    function shape = ott.shape.Shape(drag)
      % Construct a geometrical representation of the shape
      shape = ott.shape.Cylinder('radius', drag.radius, ...
          'height', drag.height, 'rotation', drag.rotation);
    end
  end

  methods % Getters/setters
    function D = get.forwardInternal(drag)

      eta = drag.viscosity;
      L = drag.height;
      R = drag.radius;
      p = L./(2*R);
      lnp = log(p);

      % Calculate coefficients for aspect ratio
      if 1/p > 0.5
        warning('ott:drag:StokesCylinder:outside_range', ...
            ['aspect ratio outside range, using aspect ratio of 0.5', ...
            newline, 'consider using StokesLmab* instead']);
        Gamma = drag.translationTable(1, 2:3);
        Delta = drag.rotationTable(1, 2:3);
      else
        Gamma = interp1(drag.translationTable(:, 1), ...
            drag.translationTable(:, 2:3), 1./p);
        Delta = interp1(drag.rotationTable(:, 1), ...
            drag.rotationTable(:, 2:3), 1./p);
      end

      % Equations from papers
      Txx = (4*eta*pi*L)./(lnp + Gamma(1));
      Tzz = (2*eta*pi*L)./(lnp + Gamma(2));
      Gxx = (pi*eta*L.^3./3)./(lnp + Delta(1));
      Gzz = (1 + Delta(2)).*(drag.A0.*pi*eta.*L.*R.^2);

      % Assemble drag tensor
      D = diag([Txx, Txx, Tzz, Gxx, Gxx, Gzz]);
    end

    function drag = set.radius(drag, val)
      assert(isnumeric(val) && isscalar(val), ...
          'radius must be numeric scalar');
      drag.radius = val;
    end

    function drag = set.height(drag, val)
      assert(isnumeric(val) && isscalar(val), ...
          'height must be numeric scalar');
      drag.height = val;
    end
  end
end

