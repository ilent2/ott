classdef Stokes < ott.utils.RotateHelper ...
    & matlab.mixin.Heterogeneous
% Base class for 6-vector force/torque drag tensors.
% Inherits from :class:`ott.utils.RotateHelper`.
%
% This class is the base class for drag tensors which can be described
% by a 3x3 translational, rotational, and cross-term matrices in
% Cartesian coordinates.
%
% Properties
%   - forward      -- Rank 2 forward tensor
%   - inverse      -- Rank 2 inverse tensor
%   - rotation     -- Rotation matrix to apply to forward/inverse
%
% Abstract properties
%   - forwardInternal   -- Forward tensor without rotation
%   - inverseInternal   -- Inverse tensor without rotation
%
% Dependent properties
%   - gamma        -- Translational component of tensor
%   - delta        -- Rotational component of tensor
%   - crossterms   -- Cross-terms component of tensor (UD)
%   - igamma       -- Inverse translational component of tensor
%   - idelta       -- Inverse rotational component of tensor
%   - icrossterms  -- Inverse cross-terms component of tensor (UD)
%
% Methods
%   - inv         -- Return the inverse or calculate the inverse
%   - mtimes      -- Multiply the tensor or calculate the forward and mul
%   - vecnorm     -- Compute vecnorm of tensor rows
%   - diag        -- Get the diagonal of the forward tensor
%   - rotate*     -- Rotate the tensor around the X,Y,Z axis
%
% Static methods
%   - FromShape   -- Construct a drag tensor from a shape or array
%
% Supported casts
%   - double            -- Get the forward drag tensor
%   - StokesData        -- Pre-computes forward/inverse tensors
%   - ott.shape.Shape  -- Gets a geometrical representation of the object
%
% See also :class:`StokesSphere` and :class:`StokesLambNn`.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    forward           % 6x6 forward drag tensor
    inverse           % 6x6 inverse drag tensor

    gamma             % Translational component of tensor
    delta             % Rotational component of tensor
    crossterms        % Cross-terms component of tensor (UD)

    igamma            % Inverse translational component of tensor
    idelta            % Inverse rotational component of tensor
    icrossterms       % Inverse cross-terms component of tensor (UD)
  end

  properties (Abstract)
    forwardInternal   % Forward drag tensor (no rotation)
    inverseInternal   % Inverse drag tensor (no rotation)
  end

  properties
    rotation          % Rotation to be applied to inverse/forward
  end

  methods (Static)
    function drag = FromShape(shape, varargin)
      % Attempt to select an appropriate drag calculation method.
      %
      % Not all shapes are supported, most shapes default to a sphere
      % whose radius matches the maxRadius of the shape.
      %
      % If the input is a single shape, attempts to choose an appropriate
      % method for the shape.  May raise warnings if no good method is found.
      %
      % If the input contains a plane and a finite extent particle,
      % models the particle as a sphere near a wall.  Raises a warning
      % if the particle is not a sphere.
      %
      % If the input contains a particle inside the circumscribing sphere
      % of another particle, uses EccentricSpheresNn.  Raises a warning
      % if the particles are not spheres.
      %
      % For all other arrays, assumes the particles are non-interacting
      % and generates drag tensors for each particle.
      % If the particles are within 5 radii of each other, raises a warning.
      %
      % Usage
      %   drag = Stokes.FromShape(shape, eta, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- A shape or array of shapes.
      %   - eta (numeric) -- Viscosity (passed to class constructor).
      %
      % Additional parameters passed to class constructors.

      if numel(shape) == 1

        if isa(shape, 'ott.shape.Sphere')
          drag = ott.drag.StokesSphere(shape.radius, varargin{:});
        elseif isa(shape, 'ott.shape.Cylinder')

          p = shape.height ./ (2*shape.radius);
          if p < 2
            drag = ott.drag.StokesStarShaped.FromShape(shape, varargin{:});
          else
            % Only works for slender cylinders
            drag = ott.drag.StokesCylinder(shape.radius, ...
                shape.height, varargin{:});
          end

        elseif shape.starShaped
          drag = ott.drag.StokesStarShaped.FromShape(shape, varargin{:});
        else

          % Choose a method based on aspect ratio
          bb = diff(shape.boundingBox);
          radius = vecnorm(bb(1:2))./2;
          height = bb(3);
          as = height./radius;

          if as > 10
            warning('ott:drag:Stokes:approx_as_cylinder', ...
              'Approximating shape as cylinder');
            drag = ott.drag.StokesCylinder(radius, ...
                height, varargin{:});
          else
            warning('ott:drag:Stokes:approx_as_sphere', ...
              'Approximating shape as sphere');
            drag = ott.drag.StokesSphere(shape.maxRadius, varargin{:});
          end
        end

      elseif numel(shape) == 2 && ( isa(shape(1), 'ott.shape.Plane') ...
          || isa(shape(2), 'ott.shape.Plane'))

        % Defer to SphereWall/FromShape
        drag = ott.drag.StokesSphereWall.FromShape(shape, varargin{:});

      elseif numel(shape) == 2 && ...
          any(vecnorm(shape(1).position - shape(2).position) ...
          < [shape.maxRadius])

        % Get inner and outer shapes
        if shape(1).maxRadius > shape(2).maxRadius
          outer = shape(1);
          inner = shape(2);
        else
          outer = shape(2);
          inner = shape(1);
        end

        if ~isa(inner, 'ott.shape.Sphere')
          warning('ott:drag:Stokes:approx_as_sphere', ...
            'Approximating particle shape as sphere');
        end

        if ~isa(outer, 'ott.shape.Sphere')
          warning('ott:drag:Stokes:approx_wall_as_sphere', ...
            'Approximating wall shape as sphere');
        end

        % Calculate separation
        separation = vecnorm(outer.position - inner.position);
        separation = outer.radius - separation - inner.radius;

        drag = ott.drag.EccentricSpheresNn(inner.maxRadius, ...
            outer.maxRadius, separation, varargin{:});
      else
        % Calculate distances between shapes
        xyz = [shape.position];
        dist = vecnorm(xyz - reshape(xyz, 3, 1, []));

        limit = 5*[shape.maxRadius];
        if any(any(dist < limit))
          warning('ott:drag:Stokes:no_interaction_terms', ...
              'Drag calculation does not include interaction terms');
        end

        % Call FromShape for each shape
        drag = ott.drag.Stokes.empty(1, 0);
        for ii = 1:numel(shape)
          drag(ii) = ott.drag.Stokes.FromShape(shape(ii), varargin{:});
        end
      end
    end
  end

  methods
    function obj = Stokes(varargin)
      % Construct a new drag tensor instance.
      %
      % Usage
      %   drag = Stokes(...)
      %
      % Optional named parameters
      %   - rotation (3x3 numeric) -- Initial rotation property.
      %     Default: ``eye(3)``.

      p = inputParser();
      p.addParameter('rotation', eye(3));
      p.parse(varargin{:});

      obj.rotation = p.Results.rotation;
    end

    function drag = ott.drag.StokesData(drag)
      % Pre-compute the forward and inverse drag tensors.
      %
      % This converts the type to a :class:`StokesData` instance using
      % the `forward` and `inverse` properties.
      %
      % If the class inherits from :class:`mixin.CalcInvDrag`, only
      % uses the `forward` property.
      %
      % Usage
      %   drag = ott.drag.StokesData(drag)

      if isa(drag, 'ott.drag.mixin.CalcInvDrag')
        drag = ott.drag.StokesData(drag.forward);
      else
        drag = ott.drag.StokesData(drag.forward, drag.inverse);
      end
    end

    function iTensor = inv(obj)
      % Return the inverse tensor if set or the inverse of the forward tensor
      %
      % Usage:
      %    inverse_tensor = inv(tensor) returns a 2-D matrix representing
      %    the tensor inverse.
      %
      % If tensor.inverse is unset, attempts to calculate the inverse from
      % the forward tensor.

      iTensor = obj.inverse;
    end

    function drag = mtimes(obj, vec)
      % Calculate the drag using the tensor and a velocity vector
      %
      % Usage:
      %   drag = tensor * vec; where vec is the velocity vector.
      %
      % If tensor.forward is not set, attempts to calculate the drag
      % from the inverse of the inverse drag tensor.

      drag = obj.forward * vec;
    end

    function num = vecnorm(obj, varargin)
      % Returns the 2-norm of each row of the forward drag tensor
      %
      % Usage
      %   N = vecnorm(obj) Calculates the vecnorm of each row
      %   of the forward drag tensor.
      %
      %   N = vecnorm(obj, p)
      %   Specifies the type of norm.  Default: 2.
      %
      %   N = vecnorm(obj, p, dim)
      %   Specifies the dimension.

      % Parse inputs
      p = inputParser;
      p.addOptional('p', 2);
      p.addOptional('dim', 2);
      p.parse(varargin{:});

      % Calculate vecnorm
      num = vecnorm(obj.forward, p.Results.p, p.Results.dim);
    end

    function mat = diag(obj)
      % Return the diagonal of the forward drag tensor
      mat = diag(obj.forward);
    end

    function mat = double(obj)
      mat = obj.forward;
    end
  end

  methods (Hidden)
    function drag = rotateInternal(drag, mat)
      % Applies a 3x3 rotation matrix to the current rotation
      % Method used by RotateHelper

      % Check for output arguments
      ott.utils.nargoutCheck(drag, nargout);

      drag.rotation = mat * drag.rotation;
    end
  end

  methods % Getters/setters
    function drag = set.rotation(drag, val)
      assert(isnumeric(val) && ismatrix(val) && all(size(val) == [3, 3]), ...
          'rotation must be 3x3 numeric matrix');
      drag.rotation = val;
    end

    function D = get.forward(drag)
      % Get drag after applying rotation

      % Get internal drag
      D = drag.forwardInternal;

      % Apply rotation to blocks
      D(1:3, 1:3) = ott.utils.rotate_3x3tensor(D(1:3, 1:3), drag.rotation.');
      D(1:3, 4:6) = ott.utils.rotate_3x3tensor(D(1:3, 4:6), drag.rotation.');
      D(4:6, 1:3) = ott.utils.rotate_3x3tensor(D(4:6, 1:3), drag.rotation.');
      D(4:6, 4:6) = ott.utils.rotate_3x3tensor(D(4:6, 4:6), drag.rotation.');
    end
    function drag = set.forward(drag, D)
      assert(ismatrix(D) && all(size(D) == [6,6]) && isnumeric(D), ...
          'inverse must be 6x6 numeric matrix');

      % Apply rotation to blocks
      D(1:3, 1:3) = ott.utils.rotate_3x3tensor(D(1:3, 1:3), drag.rotation);
      D(1:3, 4:6) = ott.utils.rotate_3x3tensor(D(1:3, 4:6), drag.rotation);
      D(4:6, 1:3) = ott.utils.rotate_3x3tensor(D(4:6, 1:3), drag.rotation);
      D(4:6, 4:6) = ott.utils.rotate_3x3tensor(D(4:6, 4:6), drag.rotation);

      % Call internal set method (may not be supported)
      drag.forwardInternal = D;
    end

    function D = get.inverse(drag)
      % Get drag after applying rotation

      % Get internal drag
      D = drag.inverseInternal;

      % Apply rotation to blocks
      D(1:3, 1:3) = ott.utils.rotate_3x3tensor(D(1:3, 1:3), drag.rotation);
      D(1:3, 4:6) = ott.utils.rotate_3x3tensor(D(1:3, 4:6), drag.rotation);
      D(4:6, 1:3) = ott.utils.rotate_3x3tensor(D(4:6, 1:3), drag.rotation);
      D(4:6, 4:6) = ott.utils.rotate_3x3tensor(D(4:6, 4:6), drag.rotation);
    end
    function drag = set.inverse(drag, D)
      assert(ismatrix(D) && all(size(D) == [6,6]) && isnumeric(D), ...
          'inverse must be 6x6 numeric matrix');

      % Apply rotation to blocks
      D(1:3, 1:3) = ott.utils.rotate_3x3tensor(D(1:3, 1:3), drag.rotation.');
      D(1:3, 4:6) = ott.utils.rotate_3x3tensor(D(1:3, 4:6), drag.rotation.');
      D(4:6, 1:3) = ott.utils.rotate_3x3tensor(D(4:6, 1:3), drag.rotation.');
      D(4:6, 4:6) = ott.utils.rotate_3x3tensor(D(4:6, 4:6), drag.rotation.');

      % Call internal set method (may not be supported)
      drag.inverseInternal = D;
    end

    function val = get.gamma(obj)
      val = obj.forward(1:3, 1:3);
    end
    function val = get.delta(obj)
      val = obj.forward(4:6, 4:6);
    end
    function val = get.crossterms(obj)
      val = obj.forward(1:3, 4:6);
    end

    function val = get.igamma(obj)
      val = obj.inverse(1:3, 1:3);
    end
    function val = get.idelta(obj)
      val = obj.inverse(4:6, 4:6);
    end
    function val = get.icrossterms(obj)
      val = obj.inverse(1:3, 4:6);
    end
  end

  methods (Static, Sealed, Access = protected)
    function default_object = getDefaultScalarElement
      % Default object is the identity tensor
      default_object = ott.drag.StokesData(eye(6), eye(6));
    end
  end
end
