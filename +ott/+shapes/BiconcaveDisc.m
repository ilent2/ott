classdef BiconcaveDisc < ott.shapes.AxisymShape & ott.shapes.Extrusion
% A biconcave disc shape.
% Inherits from :class:`Extrusion` and :class:`AxisymShape`.
%
% This shape can be used to model cells such as unstressed Red Blood Cells.
% It implements the function::
%
%   z(r) = D \sqrt{1 - \frac{4r^2}{D^2}} \left(a_0 +
%       \frac{a_1 r^2}{D^2} + \frac{a_2 r^4}{D^4} \right)
%
% where :math:`D` is the particle diameter and :math:`a` are shape
% coefficients.
%
% Methods
%   - extrudeProfile -- Calculate the height as a function of radius
%   - insideXyz      -- Determines if a point is inside, calls extrudeProfile
%   - insideRtp      -- Spherical coordinate inputs, calls insideXyz
%   - normalsRtp     -- Calculate normals at surface location
%   - normalsXyz     -- Calculate normals at surface location
%   - get_maxRadius  -- Numerically estimates maxRadius
%   - get_volume     -- Numerically estimates volume
%   - axialSymmetry  -- Particle is rotationally symmetric about z
%
% Properties
%   - radius       -- Radius of the disc
%   - coefficients -- Coefficients describing the discs shape [a0, a1, a2]
%
% Inherited properties
%   - maxRadius   -- maximum distance from shape origin
%   - maxXyRadius -- maximum radius in the xy-plane
%   - volume      -- volume of shape
%   - position    -- Location of shape ``[x, y, z]``

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    radius       % Radius of the disc
    coefficients % Coefficient describing the shape
  end

  methods (Hidden)
    function r = get_maxXyRadius(shape)
      % Get the maximum radius in the x-y plane
      r = shape.radius;
    end

    function p = get_perimiter(shape)
      % Calculate perimeter numerically

      rho = linspace(0, shape.maxRadius, 100);
      z = shape.extrudeProfile(rho);

      p = 2.0 * sum(sqrt((rho(2:end) - rho(1:end-1)).^2 ...
          + (z(2:end) - z(1:end-1)).^2));
    end
  end

  methods
    function shape = BiconcaveDisc(radius, coefficients)
      % Construct a new biconcave disc shape representation
      %
      % Usage
      %   shape = BiconcaveDisc()
      %   Constructs a disc with radius 1 and the default coefficients.
      %
      %   shape = BiconcaveDisc(radius)
      %   shape = BiconcaveDisc(radius, [a0, a1, a2])
      %
      % Parameters
      %   - radius (numeric) -- radius of the disc.
      %   - [a0, a1, a2] (numeric) -- coefficients describing the shape.
      %
      % For a shape similar to a red blood cell, use the default coefficients
      % ``a = [0.0518, 2.0026, -4.491]`` and set ``radius = 7.82``.

      % Handle default arguments
      if nargin < 2
        coefficients = [0.0518, 2.0026, -4.491];
        if nargin < 1
          radius = 1.0;
        end
      end

      shape = shape@ott.shapes.AxisymShape();

      shape.radius = radius;
      shape.coefficients = coefficients;
    end

    function z = extrudeProfile(shape, r)
      % Calculate the radial profile of the shape
      %
      % The profile for a biconcave disc is::
      %
      %   z(r) = D \sqrt{1 - \frac{4r^2}{D^2}} \left(a_0 +
      %       \frac{a_1 r^2}{D^2} + \frac{a_2 r^4}{D^4} \right)
      %
      % Where D is the shape diameter.  Values outside the disc are void.
      %
      % Usage
      %   z = shape.radialProfile(r)

      r2 = r.^2;
      a0 = shape.coefficients(1);
      a1 = shape.coefficients(2);
      a2 = shape.coefficients(3);

      z = 2.0 .* shape.radius .* sqrt(1 - r2./shape.radius.^2) ...
        .* (a0 + (a1.*r2) ./ shape.radius.^2 ./ 4 ...
        + (a2.*r2.^2) ./ shape.radius.^4 ./ 16);

      % Void values outside shape
      z(r > shape.radius) = nan;
    end

    function varargout = axialSymmetry(shape)
      % Get rotational symmetry about the primary axes.
      %
      % Particle is rotationally symmetric about z.  Returns [2, 2, 0].
      %
      % Usage
      %   [sz, sy, sz] = shape.axialSymmetry()
      %
      %   S = shape.axialSymmetry()

      s = [2, 2, 0];

      if nargout == 1
        varargout{1} = vals;
      else
        varargout{1} = vals(1);
        varargout{2} = vals(2);
        varargout{3} = vals(3);
      end
    end
  end

  methods % Getters/setters
    function shape = set.radius(shape, val)
      % Check inputs
      assert(isnumeric(val) && isscalar(val), ...
        'radius must be numeric scalar');
      shape.radius = val;
    end
    function shape = set.coefficients(shape, val)
      % Check inputs
      assert(numel(val) == 3 && isnumeric(val), ...
        'coefficients must be 3 element numeric array');
      shape.coefficients = val(:).';
    end
  end
end

