classdef Sphere < Stokes6
% Drag tensor for a sphere with Stokes Drag
%
% Properties
%   - radius    -- Radius of sphere
%   - eta       -- Viscosity of medium

  properties (SetAccess=private)
    radius        % Radius of sphere
    eta           % Viscosity of medium
  end

  methods
    function obj = Sphere(radius, eta)
      % Calculate drag tensors for spherical particle in Stokes drag.
      %
      % Usage:
      %   tensor = Sphere(radius, eta)
      %
      % Parameters
      %   - radius    -- Radius of particle
      %   - eta       -- Viscosity of medium

      obj = obj@Stokes6('translation', 6*pi*eta*radius*eye(3), ...
        'rotation', 8*pi*eta*radius.^3*eye(3), ...
        'finalize', true);

      obj.radius = radius;
      obj.eta = eta;
    end
  end
end
