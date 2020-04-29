classdef Vector
% A class describing a vector with a origin and direction
%
% Properties
%   - origin    -- Vector origins, n-dimensional array with 3 rows
%   - direction -- Vector directions, n-dimensional array with 3 rows
%
% Loosely based on Shapes.Vector from OTGO

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    origin     % Vector origins, n-dimensional array with 3 rows
    direction  % Vector directions, n-dimensional array with 3 rows
  end

  properties (Hidden, SetAccess=protected)
    data       % Actual vector data (n-dimensional with 6 rows)
  end

  methods
    function vec = Vector(origin, direction)
      % Construct a new vector instance
      %
      % Usage
      %   vec = Vector(origin, direction)
      %
      %   vec = Vector() Constructs an empty vector instance.
      %
      %   vec = Vector(direction) Constructs a vector instance with the
      %   origin of each vector set to ``[0;0;0]``.
      %
      %   vec = Vector(data) Construct a vector instance with the
      %   6 row, n-dimensional data for the origin and directions.
      %   The first three rows should describe the origin.
      %
      % Parameters
      %   - origin (numeric) -- Origins, n-dimensional array with 3 rows
      %   - direction (numeric) -- Directions, n-dimensional array with 3 rows

      % Handle no arguments
      if nargin == 0
        origin = [];
        direction = [];
      elseif nargin == 1
        if isa(origin, 'ott.utils.Vector')
          % TODO: Is there a nicer way to copy
          direction = origin.direction;
          origin = origin.origin;

        elseif size(origin, 1) == 6
          sz = size(origin);
          direction = reshape(origin(4:6, :), [3, sz(2:end)]);
          origin = reshape(origin(1:3, :), [3, sz(2:end)]);

        else
          direction = origin;
          origin = zeros(size(direction));

        end
      end

      % Store the data
      vec = vec.setData(origin, direction);
    end

    function vec = setData(vec, origin, direction)
      % Set the vector data (both origin and direction)
      %
      % Usage
      %   new_vec = vec.setData(origin, direction)
      %
      % Parameters
      %   - origin (numeric) -- Origins, n-dimensional array with 3 rows
      %   - direction (numeric) -- Directions, n-dimensional array with 3 rows

      assert(size(origin, 1) == 3, 'origin must have 3 rows');
      assert(size(direction, 1) == 3, 'direction must have 3 rows');
      assert(all(size(origin) == size(direction)), ...
          'origin and direction must have same size');
      assert(isnumeric(origin) && isreal(origin), ...
          'origin must be real numeric matrix');
      assert(isnumeric(direction) && isreal(direction), ...
          'origin must be real numeric matrix');

      vec.data = [origin; direction];
    end

    function varargout = plot(vec, varargin)
      % Plots the vector set in 3-D.
      %
      % Uses the quiver function to generate a visualisation of the
      % vector set.
      %
      % Usage
      %   h = plot(vec, ...)
      %
      % Optional named arguments
      %   - Scale (numeric) -- rescales the coordinates and components
      %     of the vector before plotting.  Can either be a scalar
      %     or vector ``[S1, S2]`` specifying separate scaling for the
      %     coordinates and components.  Default: ``[1, 1]``.
      %
      % Any unmatched named arguments are applied to the plot handle
      % returned by the quiver function.

      % Parse inputs
      p = inputParser;
      p.keepUnmatched;
      p.addParameter('Scale', [1, 1]);
      p.parse(varargin{:});

      S1 = p.Results.Scale(1);
      S2 = p.Results.Scale(2);

      % Generate plot
      h = quiver3(S1*vec.origin(1, :), S1*vec.origin(2, :), ...
          S1*vec.origin(3, :), S2*vec.direction(1, :), ...
          S2*vec.direction(2, :), S2*vec.direction(3, :), 0);

      % Apply unmatched arguments to plot handle
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];
      set(h, unmatched{:});

      % Assign outputs
      if nargout > 0
        varargout{1} = h;
      end
    end

    % TODO: disp method

    function sz = size(vec, dim)
      % Get the number of vectors contained in this object
      %
      % Usage
      %   sz = size(vec) Returns the size of the vectors.
      %
      %   sz = size(vec, dim) Returns the size of the specified dimension.

      sz = size(vec.data);
      sz = sz(2:end);

      if nargin == 2
        sz = sz(dim);
      elseif numel(sz) == 1
        sz = [sz, 1];
      end
    end

    function num = numel(vec)
      % Get the number of vectors in this object
      %
      % Usage
      %   num = numel(vec)

      num = prod(vec.size());
    end

    function vec = uminus(vec)
      % Unitary minus
      %
      % Inverts the components of the vector, leaves the origin unchanged.
      % This is somewhat equivalent to ``vec.direction = -vec.direction``.
      %
      % Usage
      %   vec = -vec
      %
      %   vec = uminus(vec)

      vec.direction = - vec.direction;
    end

    function vec = plus(vec1, vec2)
      % Binary addition
      %
      % Adds two :class:`Vector`s or a suitably sized matrix.
      % For two :class:`Vector`s, this is somewhat equivalent to
      % ``vec.direction = vec.direction + vec2.direction``.
      %
      % If one argument is not a :class:`Vector`, adds the non-vector
      % to the direction component of the vector.
      %
      % The new origin is set to the first vector.
      %
      % Usage
      %   new_vec = vec1 + vec2
      %
      %   new_vec = vec + other

      if isa(vec1, 'ott.utils.Vector') && isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, vec1.direction + vec2.direction);
      elseif isa(vec1, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, vec1.direction + vec2);
      elseif isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec2.origin, vec2.direction + vec1);
      else
        error('Expected at least one Vector instance');
      end
    end

    function vec = minus(vec1, vec2)
      % Binary subtraction
      %
      % See notes in :meth:`plus`.
      %
      % Usage
      %   new_vec = vec1 - vec2
      %
      %   new_vec = vec - other
      %   new_vec = other - vec

      if isa(vec1, 'ott.utils.Vector') && isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, vec1.direction - vec2.direction);
      elseif isa(vec1, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, vec1.direction - vec2);
      elseif isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec2.origin, vec1 - vec2.direction);
      else
        error('Expected at least one Vector instance');
      end
    end

    function vec = cross(vec1, vec2)
      % Vector cross-product
      %
      % Calculates the cross-product of the vector components.
      % This replaces some of the functionality of mtimes in OTGO.
      %
      % Usage
      %   new_vec = cross(vec1, vec2)
      %
      % Parameters
      %   - vec1, vec2 -- Can be :class:`Vector` or 3 row matrices
      %     which can be crossed with ``vec.direction``.

      if isa(vec1, 'ott.utils.Vector') && isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, ...
            cross(vec1.direction, vec2.direction));
      elseif isa(vec1, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, cross(vec1.direction, vec2));
      elseif isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec2.origin, cross(vec1, vec2.direction));
      else
        error('Expected at least one Vector instance');
      end
    end

    function vec = mtimes(vec1, vec2)
      % Scalar multiplication
      %
      % Implements scalar multiplication.  For cross-products see
      % :meth:`cross`.
      %
      % Usage
      %   new_vec = vec * scalar
      %   new_vec = scalar * vec

      if isscalar(vec1) && isa(vec2, 'ott.utils.Vector')
        vec = vec2;
        vec.direction = vec1 * vec.direction;
      elseif isscalar(vec2) && isa(vec1, 'ott.utils.Vector')
        vec = vec1;
        vec.direction = vec.direction * vec2;
      else
        error('One input must be scalar the other a Vector');
      end
    end

    function vec = times(vec1, vec2)
      % Array multiplication (element-by-element multiplication)
      %
      % For vector scalar product, see :meth:`dot`.
      %
      % Usage
      %   new_vec = vec * vec

      if isa(vec1, 'ott.utils.Vector') && isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, ...
            vec1.direction .* vec2.direction);
      elseif isa(vec1, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, vec1.direction .* vec2);
      elseif isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec2.origin, vec1 .* vec2.direction);
      else
        error('Expected at least one Vector instance');
      end
    end

    function m = dot(vec1, vec2)
      % Vector scalar product
      %
      % Usage
      %   new_vec = vec1.dot(vec2)
      %
      %   new_vec = dot(vec1, vec2)

      if isa(vec1, 'ott.utils.Vector') && isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, ...
            dot(vec1.direction, vec2.direction));
      elseif isa(vec1, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, dot(vec1.direction, vec2));
      elseif isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec2.origin, dot(vec1, vec2.direction));
      else
        error('Expected at least one Vector instance');
      end
    end

    % TODO: Outer-product?

    function vec = rdivide(vec1, vec2)
      % Right array divide (element-by-element division)
      %
      % Usage
      %   new_vec = vec1 ./ vec2
      %
      %   new_vec = vec ./ other
      %   new_vec = other ./ vec

      if isa(vec1, 'ott.utils.Vector') && isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, ...
            vec1.direction ./ vec2.direction);
      elseif isa(vec1, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec1.origin, vec1.direction ./ vec2);
      elseif isa(vec2, 'ott.utils.Vector')
        vec = ott.utils.Vector(vec2.origin, vec1 ./ vec2.direction);
      else
        error('Expected at least one Vector instance');
      end
    end

    % TODO: Normalisation

  end

  methods
    function vec = set.origin(vec, val)

      assert(size(val, 1) == 3, 'origin must have 3 rows');
      assert(all(size(val) == size(vec.direction)), ...
          'origin and direction must have same size');
      assert(isnumeric(val) && isreal(val), ...
          'origin must be real numeric matrix');

      vec.data(1:3, :) = val(:, :);
    end
    function val = get.origin(vec)
      val = reshape(vec.data(1:3, :), [3, vec.size()]);
    end

    function vec = set.direction(vec, val)

      assert(size(val, 1) == 3, 'direction must have 3 rows');
      assert(all(size(vec.origin) == size(val)), ...
          'origin and direction must have same size');
      assert(isnumeric(val) && isreal(val), ...
          'direction must be real numeric matrix');

      vec.data(4:6, :) = val(:, :);
    end
    function val = get.direction(vec)
      val = reshape(vec.data(4:6, :), [3, vec.size()]);
    end
  end
end
