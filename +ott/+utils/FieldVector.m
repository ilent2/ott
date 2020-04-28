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
        loc = [];
        vec = [];
        basis = 'cartesian';
      end

      field = field@double(vec);
      field.locations = loc;
      field.basis = basis;
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
