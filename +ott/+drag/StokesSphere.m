classdef StokesSphere < ott.drag.Stokes
% Drag tensor for a sphere with Stokes Drag
%
% Properties
%   - radius    -- Radius of sphere
%   - viscosity -- Viscosity of medium

  properties (SetAccess=protected)
    radius        % Radius of sphere
  end

  methods (Static)
    function drag = simple(shape, varargin)
      % Calculate the drag using the Stokes sphere approximation.
      %
      % Usage:
      %   drag = ott.drag.Sphere.simple(shape, ...) construct a new
      %   drag tensor for the given ott.shapes.Shape object.
      %   Shape must implement a maxRadius method, the result is used
      %   as the sphere radius.
      %
      %   drag = ott.drag.Sphere.simple(name, parameters, ...) constructs
      %   a new shape described by name and parameters.
      %   See ott.shapes.Shape.simple for supported shapes.
      %   Constructed shape must have a maxRadius method.
      %
      % See ott.drag.Sphere/Sphere for optional arguments.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('parameters', []);
      p.addParameter('viscosity', 1.0);
      p.parse(varargin{:});

      % Get a shape object from the inputs
      if ischar(shape) && ~isempty(p.Results.parameters)
        shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
        varargin = varargin(2:end);
      elseif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
        error('Must input either Shape object or string and parameters');
      end

      % Get the shape radius
      radius = shape.maxRadius;

      % Calculate drag
      drag = ott.drag.StokesSpehre(radius, p.Results.viscosity);
    end
  end

  methods
    function obj = StokesSphere(radius, varargin)
      % Calculate drag tensors for spherical particle in Stokes drag.
      %
      % Usage:
      %   tensor = Sphere(radius, eta, ...)
      %
      % Parameters
      %   - radius    -- Radius of particle
      %   - viscosity -- Viscosity of medium (optional, default: 1.0)
      %
      % Optional named arguments:
      %   - finalize (logical) -- calculate inverse drag tensor.
      %     Default: `true`.
      
      p = inputParser;
      p.addOptional('viscosity', 1.0, @(x)isnumeric(x)&&isscalar(x));
      p.addParameter('finalize', true);
      p.parse(varargin{:});

      obj = obj@ott.drag.Stokes(...
        'translation', 6*pi*p.Results.viscosity*radius*eye(3), ...
        'rotation', 8*pi*p.Results.viscosity*radius.^3*eye(3), ...
        'finalize', p.Results.finalize, 'viscosity', p.Results.viscosity);

      obj.radius = radius;
    end
  end
end
