classdef FaxenSphere < ott.drag.StokesSphereWall
% Stokes drag with Faxen's corrections for movement near a plane.
% Inherits from :class:`StokesSphereWall`.
%
% Faxen's correction can provide reasonable estimates for the drag
% on a spherical particle moving near a planar surface.
% The approximation works well to within about 1 particle radius from
% the surface.
%
% This implementation uses the Faxen's corrections described in
%
%   J. Leach, et al. Phys. Rev. E 79, 026301
%   https://doi.org/10.1103/PhysRevE.79.026301
%
% For the rotational-translational Faxen's coupling, we use Eq. 7-4.29 from
%
%   John Happel and Howard Brenner,
%   Low Reynolds number hydrodynamics, (1983)
%   https://doi.org/10.1007/978-94-009-8352-6
%
% This class assumes the surface is perpendicular to the z axis, positioned
% bellow the spherical particle.
%
% Properties
%   - radius      -- Radius of the sphere
%   - viscosity   -- Viscosity of the medium
%   - separation  -- Distance between the sphere centre and plane
%   - forward     -- Computed drag tensor
%   - inverse     -- Computed from `forward`.
%
% See :class:`Stokes` for other methods/parameters.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    forwardInternal
  end

  methods
    function drag = FaxenSphere(varargin)
      % Construct a new Faxen's corrected sphere drag tensor.
      %
      % Usage
      %   drag = FaxenSphere(radius, separation, viscosity, ...)
      %   Radius and separation should be specified in the same units.
      %
      % Parameters
      %   - radius     -- (numeric) Radius of sphere (default: 1)
      %   - separation -- Separation between sphere centre and surface
      %   - viscosity -- Viscosity of medium (default: 1.0)
      %
      % Parameters can also be passed as named arguments.
      % Unmatched parameters are passed to :class:`Stokes`.

      % Only need a constructor for help/doc functionality
      drag = drag@ott.drag.StokesSphereWall(varargin{:});
    end
  end

  methods % Getters/setters
    function D = get.forwardInternal(drag)

      % Calculate stokes sphere drag
      D = ott.drag.StokesSphere(drag.radius, drag.viscosity).forward;

      as = drag.radius ./ drag.separation;

      % Check separation is large enough
      % This threshold is determined by a 10% accuracy comparing to Chaoui
      if (drag.separation / drag.radius - 1) < 0.6
        warning('ott:drag:FaxenSphere:small_epsilon', ...
          ['Apprxomation may be inacurate for small separation', ...
          newline, 'Consider using Pade or Chaoui approximation']);
      end

      % Equation 6 and 7 from Leach
      gammaP = 1./(1 - (9/16)*as + (1/8)*as^3);
      gammaS = 1./(1 - (9/8)*as + (1/2)*as^3);

      % Equation 8 and 9 from Leach
      betaP = 1./(1 - (1/8)*as^3);
      betaS = 1./(1 - (5/16)*as^3 + (15/256)*as^6);

      % Rotation-translation coupling (Eq. 7-4.29 from Happel)
      omega = (3./32) .* as.^4 .* (1 - (3./8) .* as);

      % Apply corrections
      D(1, 1) = D(1, 1)*gammaP;
      D(2, 2) = D(2, 2)*gammaP;
      D(3, 3) = D(3, 3)*gammaS;
      D(5, 1) = -8*pi*drag.viscosity*drag.radius^2*omega;
      D(4, 2) = 8*pi*drag.viscosity*drag.radius^2*omega;
      D(4, 4) = D(4, 4)*betaS;
      D(5, 5) = D(5, 5)*betaS;
      D(6, 6) = D(6, 6)*betaP;
      D(1, 5) = (3/4) .* D(5, 1);
      D(2, 4) = (3/4) .* D(4, 2);
    end
  end
end

