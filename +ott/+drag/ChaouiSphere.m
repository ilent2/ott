classdef ChaouiSphere < ott.drag.StokesSphere
% Creeping flow around a sphere in shear flow close to a wall.
%
% This class implements a polynomial fit to the exact solution
% for spherical particles moving near a planar surface.
% The approximation should work for spacing between the sphere surface
% and plane between :math:`10^{-6}\times`radius and 1 radius.
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

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    separation     % Distance between the sphere centre and plane
  end

  methods
    function obj = ChaouiSphere(radius, separation, varargin)
      % Construct a new creeping flow sphere-wall drag tensor.
      %
      % Usage:
      %   drag = Chaoui2002Sphere(radius, separation, viscosity, ...)
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

      % l/a: ratio of the distance of the sphere centre from the
      % wall to the sphere radius
      as = separation ./ radius;

      % Calculate gap width (epsilon): as - 1
      epsilon = as - 1;

      % Check epsilon range
      % This threshold is determines from the text in Chaoui
      if epsilon > 0.2
        warning('ott:drag:ChaouiSphere:large_epsilon', ...
          ['Apprxomation may be inacurate for large separation', ...
          newline, 'Consider using Faxen or Pade approximation']);
      end

      % TODO: Should we include shear flow?
      % TODO: What about translations perpendicular to the surface?

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
      obj.forward(1, 1) = obj.forward(1, 1) .* fxx;
      obj.forward(2, 2) = obj.forward(2, 2) .* fxx;
      obj.forward(4, 4) = obj.forward(4, 4) .* cyy;
      obj.forward(5, 5) = obj.forward(5, 5) .* cyy;

      % Add cross-terms (cyx)
      obj.forward(1, 5) = -6*pi*p.Results.viscosity*radius^2*cyx;
      obj.forward(2, 4) = 6*pi*p.Results.viscosity*radius^2*cyx;
      obj.forward(5, 1) = -8*pi*p.Results.viscosity*radius^2*cyx;
      obj.forward(4, 2) = 8*pi*p.Results.viscosity*radius^2*cyx;

      % czz should have almost no effect

      % If no inverse/forward, calculate them
      if p.Results.finalize
        obj = obj.finalize();
      end
    end
  end
end

