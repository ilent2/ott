classdef (Abstract) StokesSphereWall < ott.drag.Stokes ...
    & ott.drag.mixin.VarViscosity ...
    & ott.drag.mixin.CalcInvDrag
% Abstract base class for sphere-wall drag calculation methods.
% Inherits from :class:`Stokes`.
%
% Properties
%   - radius      -- Radius of sphere
%   - separation  -- Distance from sphere centre to wall
%   - viscosity   -- Viscosity of medium
%   - inverse     -- Calculated from `forward`
%
% Abstract properties
%   - forward     -- Drag tensor calculated by method
%
% Static methods
%   - FromShape   -- Construct a drag tensor from a shape array
%
% Supported casts
%   - ott.shape.Shape -- Cast to plane and sphere

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    radius
    separation
  end
  
  properties (Hidden, SetAccess=protected)
    % Initial values only for validation (overwritten on construct)
    radiusInternal = 0
    separationInternal = Inf
  end

  methods (Static)
    function drag = FromShape(shape, varargin)
      % Attempt to choose a sphere-wall drag method based on a shape array
      %
      % Usage
      %   drag = StokesSphereWall.FromShape(shape, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- Array of shapes.  Must be two
      %     elements, one of which must be a plane.
      %
      % Additional parameters are passed to corresponding class constructor.

      p = inputParser;
      p.addOptional('viscosity', 1.0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      assert(numel(shape) == 2, 'shape must be 2 elements');

      % Get the plane/particle
      if isa(shape(1), 'ott.shape.Plane')
        plane = shape(1);
        particle = shape(2);
      else
        plane = shape(2);
        particle = shape(1);
      end

      assert(isfinite(particle.maxRadius), ...
          'particle radius must be finite');

      if ~isa(particle, 'ott.shape.Sphere')
        warning('ott:drag:StokesSphereWall:approx_as_sphere', ...
          'Approximating particle shape as sphere');
      end

      rotation = plane.rotation;
      height = abs(dot(particle.position - plane.position, plane.normal));

      if (height - particle.maxRadius)/particle.maxRadius < 1e-2
        % Doesn't work well bellow 1e-2
        drag = ott.drag.ChaouiSphere(particle.maxRadius, ...
            height, p.Results.viscosity, 'rotation', rotation, unmatched{:});
      else
        % Doesn't work for very large separation
        drag = ott.drag.PadeSphere(particle.maxRadius, ...
            height, p.Results.viscosity, 'rotation', rotation, unmatched{:});
      end
    end
  end

  methods
    function drag = StokesSphereWall(varargin)
      % Construct base class for sphere-wall methods
      %
      % Usage
      %   drag = drag@ott.drag.StokesSphereWall(radius,
      %   separation, viscosity, ...)
      %
      % Parameters
      %   - radius (numeric) -- Particle radius (default: 1.0)
      %   - separation (numeric) -- Separation (default: Inf)
      %   - viscosity (numeric) -- Viscosity (default: 1.0)
      %
      % Parameters can also be passed as named arguments.
      % Unmatched parameters are passed to :class:`Stokes`.

      p = inputParser;
      p.addOptional('radius', 1.0, @isnumeric);
      p.addOptional('separation', Inf, @isnumeric);
      p.addOptional('viscosity', 1.0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      drag = drag@ott.drag.Stokes(unmatched{:});
      drag.viscosity = p.Results.viscosity;
      drag.radius = p.Results.radius;
      drag.separation = p.Results.separation;
    end

    function shape = ott.shape.Shape(drag)
      % Construct a geometrical representation of the shape

      sph = ott.shape.Sphere(drag.radius, 'rotation', drag.rotation);
      pln = ott.shape.Plane('position', [0;0;1].*drag.separation);
      shape = [sph, pln];
    end
  end

  methods % Getters/setters
    function drag = set.radius(drag, val)
      assert(isnumeric(val) && isscalar(val), ...
          'radius must be numeric scalar');
      assert(val < drag.separation, ...
          'radius must be less than separation distance');
      drag.radius = val;
    end
    function val = get.radius(drag)
      val = drag.radiusInternal;
    end

    function drag = set.separation(drag, val)
      assert(isnumeric(val) && isscalar(val), ...
          'separation must be numeric scalar');
      assert(val > drag.radius, ...
          'separation must be greater than radius');
      drag.separation = val;
    end
    function val = get.separation(drag)
      val = drag.separationInternal;
    end
  end
end

