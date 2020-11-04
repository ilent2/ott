classdef FieldVector < double
% A class encapsulating field vector data
%
% This class keeps track of the coordinate system used to represent
% field vectors and provided convenient methods for converting between
% Cartesian and Spherical coordinates.
%
% Properties
%   - locations -- Location data for field vectors (can be empty)
%   - basis     -- Coordinate system basis (cartesian or spherical)
%
% Properties (possibly computed)
%   - vxyz   -- Data in Cartesian coordinates
%   - vrtp   -- Data in Spherical coordinates
%
% Methods
%   - plus        -- Add field vectors
%   - minus       -- Subtract field vectors
%   - mtimes,times      -- Scalar multiplication
%   - mrdivide,rdivide  -- Scalar division

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Implement subsref etc (and remove double(field) in getters)
%   This would be an optimisation.  But, we shouldn't allow direct
%   access to the data.  So field.vxyz(1, 3) should be the same
%   as field(1, 3) (without the cast), but field(1, 3) should be
%   disallowed.

  properties
    locations  % Location data for field vectors (can be empty)
    basis      % Coordinate system basis (cartesian or spherical)
  end

  properties (Dependent)
    vxyz      % Data in Cartesian coordinates
    vrtp      % Data in Spherical coordinates
  end

  methods
    function field = FieldVector(loc, vec, basis)
      % Construct a new FieldVector instance
      %
      % Usage
      %   field = FieldVector(xyz, vxyz, 'cartesian')
      %   field = FieldVector(rtp, vrtp, 'spherical')

      if nargin == 0
        loc = zeros(3, 0);
        vec = zeros(3, 0);
        basis = 'cartesian';
      end

      field = field@double(vec);
      field.locations = loc;
      field.basis = basis;
    end

    function vec = plus(v1, v2)
      % Add field vectors

      if isa(v1, 'ott.utils.FieldVector') && isa(v2, 'ott.utils.FieldVector')

        assert(isequaln(v1.locations, v2.locations), ...
            'Locations must match when adding field vectors');

        if strcmpi(v1.basis, v2.basis)
          % Just add data (like normal)
          vec = ott.utils.FieldVector(v1.locations, ...
              builtin('plus', v1, v2), v1.basis);
        else
          % Convert to a common basis and add
          vec = ott.utils.FieldVector(v1.locations, ...
              v1.vxyz + v2.vxyz, 'cartesian');
        end

      elseif isscalar(v1)
        % Keep scalar addition
      	vec = ott.utils.FieldVector(v2.locations, ...
            builtin('plus', v1, v2), v2.basis);
      elseif isscalar(v2)
        % Keep scalar addition
      	vec = ott.utils.FieldVector(v1.locations, ...
            builtin('plus', v1, v2), v1.basis);
      else
        error('Unsupported addition for FieldVector');
      end

    end

    function vec = times(v1, v2)
      % Scalar multiplication of vector data

      % Defer to mtimes method
      vec = v1 * v2;
    end

    function vec = mtimes(v1, v2)
      % Scalar multiplication

      if isscalar(v1)
        % Keep scalar mtimes
      	vec = ott.utils.FieldVector(v2.locations, ...
            builtin('mtimes', v1, v2), v2.basis);
      elseif isscalar(v2)
        % Keep scalar mtimes
      	vec = ott.utils.FieldVector(v1.locations, ...
            builtin('mtimes', v1, v2), v1.basis);
      else
        error('Only scalar multiplciation supported for FieldVector');
      end

    end

    function vec = rdivide(v1, s)
      % Scalar division
      vec = v1 / s;
    end

    function vec = mrdivide(v1, s)
      % Scalar division

      assert(isscalar(s) && isnumeric(s), ...
          'second element to mrdivide must be numeric scalar');
      vec = ott.utils.FieldVector(v1.locations, ...
          builtin('mrdivide', v1, s), v1.basis);
    end
  end

  methods % Getters/setters
    function field = set.locations(field, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
        'locations must be 3xN numeric matrix');
      field.locations = val;
    end
    
    function val = get.vxyz(field)
      if strcmpi(field.basis, 'cartesian')
        val = double(field);
      elseif strcmpi(field.basis, 'spherical')
        if isempty(field.locations)
          error('Need locations to apply coordinate transformation');
        end
        val = ott.utils.rtpv2xyzv(field.vrtp, field.locations);
      else
        error('Unable to convert between these bases');
      end
    end

    function val = get.vrtp(field)
      if strcmpi(field.basis, 'cartesian')
        if isempty(field.locations)
          error('Need locations to apply coordinate transformation');
        end
        val = ott.utils.xyzv2rtpv(field.vxyz, field.locations);
      elseif strcmpi(field.basis, 'spherical')
        val = double(field);
      else
        error('Unable to convert between these bases');
      end
    end

  end
end
