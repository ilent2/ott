classdef FaxenSphere < ott.drag.StokesSphere
% Stokes drag with Faxen's corrections for movement near a plane.
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

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    separation     % Distance between the sphere centre and plane
  end

  methods
    function obj = FaxenSphere(radius, separation, varargin)
      % Construct a new Faxen's corrected sphere drag tensor.
      %
      % Usage:
      %   drag = FaxenSphere(radius, separation, viscosity, ...)
      %   radius and separation should be specified in the same units.
      %
      % Parameters:
      %   - radius     -- Radius of sphere
      %   - separation -- Separation between sphere centre and surface
      %   - viscosity -- Viscosity of medium (optional, default: 1.0)
      %
      % Optional named arguments:
      %   - finalize (logical) -- calculate inverse drag tensor.
      %     Default: `true`.

      p = inputParser;
      p.addOptional('viscosity', 1.0, @(x)isnumeric(x)&&isscalar(x));
      p.addParameter('finalize', true);
      p.parse(varargin{:});

      % Construct regular sphere instance
      obj = obj@ott.drag.StokesSphere(radius, p.Results.viscosity, ...
          'finalize', false);
      obj.separation = separation;

      as = radius ./ separation;

      % Check separation is large enough
      % This threshold is determined by a 10% accuracy comparing to Chaoui
      if (separation ./ radius - 1) < 0.6
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
      obj.forward(1, 1) = obj.forward(1, 1)*gammaP;
      obj.forward(2, 2) = obj.forward(2, 2)*gammaP;
      obj.forward(3, 3) = obj.forward(3, 3)*gammaS;
      obj.forward(5, 1) = -8*pi*p.Results.viscosity*radius^2*omega;
      obj.forward(4, 2) = 8*pi*p.Results.viscosity*radius^2*omega;
      obj.forward(4, 4) = obj.forward(4, 4)*betaS;
      obj.forward(5, 5) = obj.forward(5, 5)*betaS;
      obj.forward(6, 6) = obj.forward(6, 6)*betaP;
      obj.forward(1, 5) = (3/4) .* obj.forward(5, 1);
      obj.forward(2, 4) = (3/4) .* obj.forward(4, 2);

      % If no inverse/forward, calculate them
      if p.Results.finalize
        obj = obj.finalize();
      end
    end
  end
end

