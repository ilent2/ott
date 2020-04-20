classdef EccentricSpheresNn < ott.drag.Stokes
% Calculate drag on an eccentric sphere using Gibson's NN approach.
%
% Uses the NN from
%
%   Lachlan J. Gibson, et al. Phys. Rev. E 99, 043304
%   https://doi.org/10.1103/PhysRevE.99.043304
%
% Properties
%   - innerRadius      -- Radius of inner sphere
%   - outerRadius      -- Radius of outer sphere
%   - separation       -- Minimum separation between spheres
%   - viscosity        -- Viscosity of surrounding fluid (default: 1.0)

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    innerRadius       % Radius of inner sphere
    outerRadius       % Radius of outer sphere
    separation        % Minimum separation between spheres
  end

  methods
    function obj = EccentricSpheresNn(innerRadius, outerRadius, ...
        separation, varargin)
      % Calculate the drag on an eccentric sphere.
      %
      % Usage:
      %   drag = EccentricSpheresNn(innerRadius, outerRadius, ...
      %     separation, viscosity)
      %
      % Parameters:
      %   - innerRadius -- Radius of inner sphere
      %   - outerRadius -- Radius of outer sphere
      %   - separation -- Minimum separation between inner and outer sphere
      %   - viscosity -- Viscosity of medium (optional, default: 1.0)
      %
      % Optional named arguments:
      %   - finalize (logical) -- calculate inverse drag tensor.
      %     Default: `true`.

      p = inputParser;
      p.addOptional('viscosity', 1.0, @(x)isnumeric(x)&&isscalar(x));
      p.addParameter('finalize', true);
      p.parse(varargin{:});

      % Construct regular Stokes instance
      obj = obj@ott.drag.Stokes(...
        'translation', 6*pi*p.Results.viscosity*innerRadius*eye(3), ...
        'rotation', 8*pi*p.Results.viscosity*innerRadius.^3*eye(3), ...
        'viscosity', p.Results.viscosity, 'finalize', false);
      
      % Check parameter ranges
      assert(innerRadius < outerRadius, ...
        'InnerRadius must be less than outerRadius');
      assert(separation <= outerRadius - innerRadius, ...
        'Separation must be less than outerRadius - innerRadius');
      
      % Set properties
      obj.separation = separation;
      obj.innerRadius = innerRadius;
      obj.outerRadius = outerRadius;

      % Calculate parameters for Lachlan's script
      D = separation ./ innerRadius;
      lam = innerRadius ./ outerRadius;

      [gy, fx, fxc, fz, gz] = WEES(D, lam);

      obj.forward(1, 1) = obj.forward(1, 1) .* fx;
      obj.forward(2, 2) = obj.forward(2, 2) .* fx;
      obj.forward(3, 3) = obj.forward(3, 3) .* fz;
      obj.forward(4, 4) = obj.forward(4, 4) .* gy;
      obj.forward(5, 5) = obj.forward(5, 5) .* gy;
      obj.forward(6, 6) = obj.forward(6, 6) .* gz;

      % Add cross-terms
      obj.forward(1, 5) = 6*pi*p.Results.viscosity*innerRadius^2*fxc;
      obj.forward(2, 4) = -6*pi*p.Results.viscosity*innerRadius^2*fxc;
      obj.forward(5, 1) = 8*pi*p.Results.viscosity*innerRadius^2*fxc;
      obj.forward(4, 2) = -8*pi*p.Results.viscosity*innerRadius^2*fxc;

      % If no inverse/forward, calculate them
      if p.Results.finalize
        obj = obj.finalize();
      end
    end
  end
end

