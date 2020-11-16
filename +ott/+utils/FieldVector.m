classdef FieldVector < double
% Base class for classes encapsulating field vector data.
%
% .. warning:: This class may not behave as expected in Matlab R2018a
%     and earlier.  Use with caution.
%
% .. warning:: Direct element access is not recommended with this class.
%
% Methods
%   - plus, minus, uminus, times, mtimes, rdivide, mrdivide
%   - sum     -- Add field vector components
%   - cross, dot -- Implementations of cross and dot products
%   - vecnorm -- Calculate vector magnitude
%   - vxyz    -- Data in Cartesian coordinates
%   - vrtp    -- Data in Spherical coordinates

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: This class is still a little dodgy, we really need to implement
% a subsasgn/subsref method to deal with this properly.

  methods (Access=protected)
    function field = FieldVector(vec, pos)
      % Construct a new FieldVector instance
      %
      % Usage
      %   field = FieldVector(vec, pos)
      %
      %   field = FieldVector([vec; pos])
      %
      % Parameters
      %   - vec (3xN numeric) -- Field vectors
      %
      %   - pos (3xN numeric) -- Coordinate locations (optional).
      %     If omitted, these are not stored, defaults to [0;0;0] when used.

      if nargin == 2
        data = [vec; pos];
      else
        data = vec;
      end

      assert(any(size(data, 1) == [3, 6]), ...
        'First dimension of data [vec; pos] must be 3 or 6');
      field = field@double(data);
    end
  end

  methods
    function fv = cast(fv, type, varargin)
      % Cast FieldVector to a differnt ata type/class.
      %
      % This is useful in R2018a where implicit cast dont always seem
      % to work correctly (the problem seems to be fixed in R2020b).
      %
      % Usage
      %   fv = cast(fv, 'like', other) -- Casts field vector type to
      %   match other.  If other is not a field vector, defers to
      %   builtin cast method.
      %
      %   fv = cast(fv, typename) -- Casts to specified field vector type.
      %   If typename is not a field vector type, defers to builtin method.
      
      atype = type;
      if nargin == 3
        atype = class(varargin{1});
      end
        
      % Defer to builtin for sparsity/complexity/everything else
      fv(:, :) = builtin('cast', fv, type, varargin{:});
      
      % Do field vector casts (if needed)
      if strcmpi(atype, 'ott.utils.FieldVectorCart')
        fv = ott.utils.FieldVectorCart(fv);
      elseif strcmpi(atype, 'ott.utils.FieldVectorSph')
        fv = ott.utils.FieldVectorSph(fv);
      end
    end
    
    function fv = vxyz(fv, varargin)
      % Get Cartesian field vector instance

      if ~isa(fv, 'ott.utils.FieldVectorCart')
        fv = ott.utils.FieldVectorCart(fv);
      end

      % Get only vector component
      sz = size(fv);
      fv = reshape(double(fv(1:3, :)), [3, sz(2:end)]);

      % Interpret additional arguments as subscripts
      if nargin >= 2
        fv = fv(varargin{:});
      end
    end

    function fv = vrtp(fv, varargin)
      % Get Spherical field vector instance

      if ~isa(fv, 'ott.utils.FieldVectorSph')
        fv = ott.utils.FieldVectorSph(fv);
      end

      % Get only vector component
      sz = size(fv);
      fv = reshape(double(fv(1:3, :)), [3, sz(2:end)]);

      % Interpret additional arguments as subscripts
      if nargin >= 2
        fv = fv(varargin{:});
      end
    end

    function out = sum(vec, dim, fvkeep)
      % Sum field vectors along specified dimension
      %
      % Usage
      %   fv = sum(fv, dim) -- Returns a field vector without any
      %   coordinate information.
      %
      %   fv = sum(fv, dim, keep) -- When ``keep`` is 'keepFirst',
      %   keeps the coordinates of the first element in the specified
      %   dimension.
      %
      % If dimension is 1, returns a double.  Otherwise returns a
      % field vector.
      
      % Cast to Cartesian
      vec = ott.utils.FieldVectorCart(vec);

      % Do sum on double data
      out = sum(vec.vxyz, dim);

      if size(out, 1) == 3
        
        % Set coordinates (if requested)
        if nargin == 3
          assert(strcmpi(fvkeep, 'keepFirst'), ...
            'keep argument must be ''keepFirst''');
          idx = repmat({':'}, 1, ndims(vec));
          idx{dim} = 1;
          out = [out; double(vec(4:6, idx{2:end}))];
        end
        
        out = ott.utils.FieldVectorCart(out);
      end
    end
    
    function vec = cross(v1, v2)
      % Vector cross product
      %
      % Usage
      %   fv = cross(fv, fv)
      
      assert(isa(v1, 'ott.utils.FieldVector') ...
        && isa(v2, 'ott.utils.FieldVector'), ...
        'Both inputs must be field vectors');
      
      vec = cross(v1.vxyz, v2.vxyz);
      vec = ott.utils.FieldVectorCart(vec);
    end
    
    function vec = dot(v1, v2)
      % Vector dot product
      %
      % Usage
      %   s = dot(fv, fv)
      
      assert(isa(v1, 'ott.utils.FieldVector') ...
        && isa(v2, 'ott.utils.FieldVector'), ...
        'Both inputs must be field vectors');
      
      vec = dot(v1.vxyz, v2.vxyz);
    end
    
    function vec = vecnorm(vec)
      % Calculate vector magnitude
      %
      % Usage
      %   s = vecnorm(fv)
      
      vec = vecnorm(vec.vxyz);
    end
    
    function vec = conj(vec)
      % Complex conjugate
      %
      % Usage
      %   fv = conj(fv)
      
      vec(1:3, :) = conj(double(vec(1:3, :)));
    end

    function vec = plus(v1, v2)
      % Addition of field vectors
      %
      % Usage
      %   fv = fv1 + fv2 -- Adds field vectors in Cartesian basis.
      %   Resulting field vector object has no location data.
      %
      %   s = fv1 + s2, and s = s1 + fv1 -- Add non-field vector types.
      %   Results in a double or other non-field vector type.

      if isa(v1, 'ott.utils.FieldVector') && isa(v2, 'ott.utils.FieldVector')
        vec = ott.utils.FieldVectorCart(v1.vxyz + v2.vxyz);
      elseif isa(v1, 'ott.utils.FieldVector')
        vec = v1.vxyz + v2;
      else
        vec = v1 + v2.vxyz;
      end
    end

    function vec = minus(v1, v2)
      % Subtraction of field vectors
      %
      % Usage
      %   fv = fv1 - fv2 -- Adds field vectors in Cartesian basis.
      %   Resulting field vector object has no location data.
      %
      %   s = fv1 - s2, and s = s1 - fv1 -- Add non-field vector types.
      %   Results in a double or other non-field vector type.

      if isa(v1, 'ott.utils.FieldVector') && isa(v2, 'ott.utils.FieldVector')
        vec = ott.utils.FieldVectorCart(v1.vxyz - v2.vxyz);
      elseif isa(v1, 'ott.utils.FieldVector')
        vec = v1.vxyz - v2;
      else
        vec = v1 - v2.vxyz;
      end
    end

    function vec = uminus(vec)
      % Unitary minus of field vector
      %
      % Usage
      %   vec = -vec -- Negates vector components.  Leaves location unchanged.

      sz = size(vec);
      vec(1:3, :) = -double(vec(1:3, :));
      vec = reshape(vec, sz);
    end

    function vec = times(v1, v2)
      % Multiplication of field vectors
      %
      % Usage
      %   fv = fv1 .* fv2 -- Times field vectors in Cartesian basis.
      %   Resulting field vector object has no location data.
      %
      %   fv = fv1 .* s2, and fv = s1 .* fv1 -- Scale field vector.
      %   Resulting field vector has same location as original.

      if isa(v1, 'ott.utils.FieldVector') && isa(v2, 'ott.utils.FieldVector')
        vec = ott.utils.FieldVectorCart(v1.vxyz .* v2.vxyz);
      elseif isa(v1, 'ott.utils.FieldVector')
        sz = size(v1);
        v1(1:3, :) = double(v1(1:3, :)) .* v2;
        vec = reshape(v1, sz);
      else
        sz = size(v2);
        v2(1:3, :) = v1 .* double(v2(1:3, :));
        vec = reshape(v2, sz);
      end
    end

    function vec = mtimes(v1, v2)
      % Matrix and scalar multiplication
      %
      % Usage
      %   s = M * fv -- Matrix multiplication.  M should be nx3.
      %   First gets the Cartesian field vectors.
      %
      %   fv = fv1 * s2, and fv = s1 * fv1 -- Scale field vector.
      %   Same behaviour as .* operation.

      if ~isa(v1, 'ott.utils.FieldVector') && ~isscalar(v1) ...
          && isa(v2, 'ott.utils.FieldVector')
        vec = v1 * v2.vxyz;
      else
        vec = v1 .* v2;
      end
    end

    function vec = rdivide(v1, s)
      % Scalar division
      %
      % Usage
      %   fv = fv ./ s -- Divide field vector values by scalar.

      assert(isa(v1, 'ott.utils.FieldVector'), ...
          'First argument must be a field vector');

      sz = size(v1);
      v1(1:3, :) = double(v1(1:3, :)) ./ s;
      vec = reshape(v1, sz);
    end

    function vec = mrdivide(v1, s)
      % Scalar division
      %
      % Usage
      %   fv = fv / s -- uses same behaviour as ./ operator.

      vec = v1 ./ s;
    end
  end
end
