classdef FaxenSphere < ott.drag.StokesSphere
% Stokes drag with Faxen's corrections for movement near a plane.
%
% Faxen's correction can provide reasonable estimates for the drag
% on a spherical particle moving near a planar surface.
% The approximation works well to within about 1 particle radius from
% the surface.
%
% This implementation uses the Faxen's corrections described in:
%
%   J. Leach, et al. Phys. Rev. E 79, 026301
%   https://doi.org/10.1103/PhysRevE.79.026301
%
% This class assumes the surface is perpendicular to the z axis, positioned
% bellow the spherical particle.
%
% Properties:
%   - radius      -- Radius of the sphere
%   - viscosity   -- Viscosity of the medium
%   - separation  -- Distance between the sphere surface and plane

  properties (SetAccess=protected)
    separation     % Distance between the sphere surface and plane
  end

  methods (Static)
    function drag = simple()
    end
  end

  methods
    function obj = FaxenSphere(radius, viscosity, separation)

      % Construct regular sphere instance
      obj = obj@StokesSphere(radius, viscosity);

      as = radius ./ separation;

      gammaP = 1./(1 - (9/16)*as + (1/8)*as^3);
      gammaS = 1./(1 - (9/8)*as + (1/2)*as^3);

      betaP = 1./(1 - (1/8)*as^3);
      betaS = 1./(1 - (5/16)*as^3 + (15/256)*as^6);

      % Apply corrections
      obj.forward(1, 1) = obj.forward(1, 1)*gammaP;
      obj.forward(2, 2) = obj.forward(2, 2)*gammaP;
      obj.forward(3, 3) = obj.forward(3, 3)*gammaS;
      obj.forward(4, 4) = obj.forward(4, 4)*betaS;
      obj.forward(5, 5) = obj.forward(5, 5)*betaS;
      obj.forward(6, 6) = obj.forward(6, 6)*betaP;

      % Compute inverse
      obj.inverse = inv(obj.forward);

    end
  end
end

