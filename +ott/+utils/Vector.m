classdef Vector < ott.utils.RotateHelper
% A class describing a vector with a origin and direction
% Inherits from :class:`RotateHelper`.
%
% Properties
%   - origin    -- Vector origins, n-dimensional array with 3 rows
%   - direction -- Vector directions, n-dimensional array with 3 rows
%
% Methods
%   - cat       -- Concatenate two vector objects
%   - rotate    -- Rotate the vectors (and optionally, their origins)
%   - rotate*   -- Rotate the vectors around the X,Y,Z axis
%   - isnumeric -- Returns true if class is Vector
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
    function vec = Vector(varargin)
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
      %   - origin (3x... numeric) -- Vector origins
      %   - direction (3x... numeric) -- Vector directions
      
      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('direction', [], @isnumeric);
      p.parse(varargin{:});
      
      if size(p.Results.origin, 1) == 6 && isempty(p.Results.direction)
        origin = p.Results.origin(1:3, :);
        direction = p.Results.origin(5:6, :);
      else
        if any(strcmpi(p.UsingDefaults, 'origin'))
          if isempty(p.Results.direction)
            origin = [];
          else
            origin = [0;0;0];
          end
        else
          origin = p.Results.origin;
        end
        direction = p.Results.direction;
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
      
      assert((ndims(origin) == ndims(direction) ...
        && all(size(origin) == size(direction))) ...
        || numel(origin) == 3 || numel(direction) == 3, ...
        'origin and direction must have same size or length 1');
      
      assert(isnumeric(origin) && isreal(origin), ...
          'origin must be real numeric matrix');
      assert(isnumeric(direction) && isreal(direction), ...
          'origin must be real numeric matrix');
      
      % Duplicate elements as required
      if numel(origin) == 3 && numel(direction) ~= 3
        sz = size(direction);
        origin = repmat(origin, [1, sz(2:end)]);
      elseif numel(direction) == 3 && numel(origin) ~= 3
        sz = size(origin);
        direction = repmat(direction, [1, sz(2:end)]);
      end

      vec.data = [double(origin); double(direction)];
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
      p.KeepUnmatched = true;
      p.addParameter('Scale', [1, 1]);
      p.parse(varargin{:});

      S1 = p.Results.Scale(1);
      S2 = p.Results.Scale(2);

      % Generate plot
      h = quiver3(S1*vec.origin(1, :), S1*vec.origin(2, :), ...
          S1*vec.origin(3, :), S2*vec.direction(1, :), ...
          S2*vec.direction(2, :), S2*vec.direction(3, :), 0);

      % Apply unmatched arguments to plot handle
      unmatched = ott.utils.unmatchedArgs(p);
      if ~isempty(unmatched)
        set(h, unmatched{:});
      end

      % Assign outputs
      if nargout > 0
        varargout{1} = h;
      end
    end

    function b = isnumeric(vec)
      % Returns true (only if class is ott.utils.Vector
      %
      % If you derive a class from Vector and the type should be
      % numeric, overload this function.

      b = strcmpi(class(vec), 'ott.utils.Vector');
    end
    
    function direction = double(vec)
      % Cast to double discards origin information
      direction = vec.direction;
    end

    function varargout = size(vec, varargin)
      % Get the number of vectors contained in this object.
      %
      % The leading dimension size is always 3 to make the class
      % more compatible as a drop-in replacement for 3xN arrays.
      %
      % Usage
      %   sz = size(vec) Returns the size of the vectors in a size array.
      %
      %   sz = size(vec, dim) Returns the size of the specified dimension.
      %
      % Supports other arguments that can be passed to builtin size.

      sz = size(vec.data);
      sz(1) = 3;
      
      [varargout{1:nargout}] = ott.utils.size_helper(sz, varargin{:});
    end

    function n = ndims(vec)
      % Calculate the number of dimensions
      %
      % Usage
      %   n = ndims(vec)   or   n = vec.ndims()

      n = length(size(vec));
    end

    function out = sum(vec, dim)
      % Calculate the sum of the vector along the specified dimension
      %
      % Usage
      %   d = sum(vec) Calculates the sum along the first dimension.
      %   Returns a double.
      %
      %   d = sum(vec, dim) Calculates the sum along the specified
      %   dimension.  Return a double.

      if nargin == 1
        dim = 1;
      end

      out = sum(vec.direction, dim);
    end
    
    function out = vecnorm(vec, type, dim)
      % Calculate the vecnorm of the vector
      %
      % Usage
      %   out = vecnorm(vec) Calculates the 2-norm along the first
      %   dimension.  Returns a double.
      %
      %   out = vecnorm(vec, type, dim) Specifies the dimension to
      %   operate over and the type.  Returns a double.
      
      if nargin < 3
        dim = 1;
        if nargin < 2
          type = 2;
        end
      end
      
      out = vecnorm(vec.direction, type, dim);
    end

    function vec = cat(dim, vec1, varargin)
      % Concatenates two or more Vector objects
      %
      % Usage
      %   vec = cat(dim, vec1, vec2, vec3, ...)
      %
      % Parameters
      %   - dim (numeric) -- Dimension to concatenate along, starting from
      %     1 for horizontal concatenation.

      vec = vec1;
      for ii = 1:length(varargin)
        ovec = varargin{ii};
        vec.data = cat(dim+1, vec.data, ovec.data);
      end
    end

    function vec = horzcat(varargin)
      % Horizontal concatenation of vectors
      %
      % Usage
      %   vec = [vec1, vec2, vec3];

      vec = cat(1, varargin{:});
    end

    function vertcat(varargin)
      % Disallow vertcat on Vectors (ambiguous)
      error('Vector only supports hozcat, use cat for other dimensions');
    end

    function num = numel(vec)
      % Get the number of vectors in this object
      %
      % Usage
      %   num = numel(vec)

      sz = vec.size();
      num = prod(sz(2:end));
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
      % Scalar and matrix multiplication
      %
      % Implements scalar multiplication.  For cross-products see
      % :meth:`cross`.
      %
      % Usage
      %   new_vec = vec * scalar
      %   new_vec = scalar * vec
      %   Scalar multiplication.
      %
      %   new_vec = R * vec
      %   Matrix-vector multiplication.  R must be a 3xN matrix.
      %   Rotates vectors about their origins, for rotation of origins
      %   see :meth:`rotate`.

      if isscalar(vec1) && isa(vec2, 'ott.utils.Vector')
        vec = vec2;
        vec.direction = vec1 * vec.direction;
      elseif isscalar(vec2) && isa(vec1, 'ott.utils.Vector')
        vec = vec1;
        vec.direction = vec.direction * vec2;
      elseif ismatrix(vec1) && all(size(vec1) == [3, 3]) ...
          && isa(vec2, 'ott.utils.Vector')
        vec = vec2;
        vec.direction = vec1 * vec.direction;
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
        m = dot(vec1.direction, vec2.direction);
      elseif isa(vec1, 'ott.utils.Vector')
        m = dot(vec1.direction, vec2);
      elseif isa(vec2, 'ott.utils.Vector')
        m = dot(vec1, vec2.direction);
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

  methods (Hidden)
    function vec = rotateInternal(vec, R, varargin)
      % Applies a rotation to vectors or vectors and origins
      %
      % Usage
      %   vec = vec.rotate(vec, R)
      %   Rotates vectors about origins.
      %
      %   vec = vec.rotate(vec, R, 'origin', true)
      %   Rotate vectors and origins.

      p = inputParser;
      p.addParameter('origin', false);
      p.parse(varargin{:});

      % Apply rotation to vectors
      vec = R * vec;

      % Apply rotation to origin too
      if p.Results.origin
        vec.origin = R * vec.origin;
      end
    end
  end

  methods % Getters/setters
    function vec = set.origin(vec, val)

      assert(size(val, 1) == 3, 'origin must have 3 rows');
      assert(all(size(val) == size(vec.direction)), ...
          'origin and direction must have same size');
      assert(isnumeric(val) && isreal(val), ...
          'origin must be real numeric matrix');

      vec.data(1:3, :) = val(:, :);
    end
    function val = get.origin(vec)
      sz = size(vec.data);
      val = reshape(vec.data(1:3, :), [3, sz(2:end)]);
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
      sz = size(vec.data);
      val = reshape(vec.data(4:6, :), [3, sz(2:end)]);
    end
  end
end
