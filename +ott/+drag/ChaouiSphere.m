classdef ChaouiSphere < ott.drag.StokesSphereWall
% Creeping flow around a sphere close to a wall.
% Inherits from :class:`StokesSphereWall`.
%
% This class implements a polynomial fit to the exact solution
% for spherical particles moving in creeping flow near a planar surface.
% The approximation should work for spacing between the sphere surface
% and plane between :math:`10^{-6}` radius and 1 radius.
%
% This class implements the polynomial approximation described in
%
%   M. Chaoui and F. Feuillebois,
%   Creeping Flow Around a Sphere in a Shear Flow Close to a Wall.
%   The Quarterly Journal of Mechanics and Applied Mathematics,
%   Volume 56, Issue 3, August 2003, Pages 381--410,
%   https://doi.org/10.1093/qjmam/56.3.381
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
    function drag = ChaouiSphere(varargin)
      % Construct a new creeping flow sphere-wall drag tensor.
      %
      % Usage:
      %   drag = ChaouiSphere(radius, separation, viscosity, ...)
      %   radius and separation should be specified in the same units.
      %
      % Parameters:
      %   - radius     -- Radius of sphere
      %   - separation -- Separation between sphere centre and surface
      %   - viscosity -- Viscosity of medium (optional, default: 1.0)
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

      % l/a: ratio of the distance of the sphere centre from the
      % wall to the sphere radius
      as = drag.separation ./ drag.radius;

      % Calculate gap width (epsilon): as - 1
      epsilon = as - 1;

      % Check epsilon range
      % This threshold is determines from the text in Chaoui
      if epsilon > 0.2
        warning('ott:drag:ChaouiSphere:large_epsilon', ...
          ['Apprxomation may be inacurate for large separation', ...
          newline, 'Consider using Faxen or Pade approximation']);
      end

      % Values from table 11
      phis = [-8/15, -64/375, 0.011595712294862, -0.002559314461340, ...
          0.002165777707452, 0.000351260314552];
      fs = [0.954293724783876, 0.429450132564500, -0.001897844702304, ...
          0.002058408405495, 0.000096108639584, -0.001248147281379];
      gammas = [-1/10, -43/250, -0.036913066460225, 0.001486892317125, ...
          0.000012689734456, 0.000103798994187];
      cs = [-0.192952745666190, 0.100579155700110, 0.094493729126963, ...
          0.003821112414990, -0.000819028830091, -0.000097511506358];

      % Equation 5.4 and 5.5
      fxx = sum(phis.*epsilon.^(0:5).*log(epsilon) + fs.*epsilon.^(0:5));
      cyx = sum(gammas.*epsilon.^(0:5).*log(epsilon) + cs.*epsilon.^(0:5));

      % Values from table 13
      gammar = [-2/5, -0.528001276176667, -0.212879560114862, ...
          -0.035965644690736, -0.006385459746252, 0.000167620439255];
      cr = [0.370892565890165, 0.340079923061464, 0.225531274283815, ...
          0.097897336215370, 0.005878696055717, 0.001503759398496];

      % Equation 5.6
      cyy = sum(gammar.*epsilon.^(0:5).*log(epsilon) + cr.*epsilon.^(0:5));

      % Modify diagonal terms (cyy, fxx)
      D(1, 1) = D(1, 1) .* fxx;
      D(2, 2) = D(2, 2) .* fxx;
      D(4, 4) = D(4, 4) .* cyy;
      D(5, 5) = D(5, 5) .* cyy;

      % Add cross-terms (cyx)
      D(1, 5) = -6*pi*drag.viscosity*drag.radius^2*cyx;
      D(2, 4) = 6*pi*drag.viscosity*drag.radius^2*cyx;
      D(5, 1) = -8*pi*drag.viscosity*drag.radius^2*cyx;
      D(4, 2) = 8*pi*drag.viscosity*drag.radius^2*cyx;

      % Add warning about using Faxen for D(3, 3) and D(6, 6)
      warning('ott:drag:ChaouiSphere:faxen_perp_terms', ...
          'Using Faxen corrections for D(3, 3) and D(6, 6)');
      gammaS = 1./(1 - (9/8)./as + (1/2)*as^(-3));
      betaP = 1./(1 - (1/8)*as^(-3));
      D(3, 3) = D(3, 3)*gammaS;
      D(6, 6) = D(6, 6)*betaP;
    end
  end
end

