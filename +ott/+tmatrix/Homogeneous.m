classdef Homogeneous < ott.tmatrix.Tmatrix
% Base class for homogeneous T-matrices
%
% Properties
%   - index_relative    -- Relative refractive index of homogeneous particle
%
% See base class for additional methods/properties.

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    index_relative      % Relative refractive index of homogeneous particle
  end

  methods
    function tm = Homogeneous(varargin)
      % Construct a new Homogeneous T-matrix object.
      %
      % Usage
      %   tmatrix = Tmatrix(...)
      %   New empty T-matrix.  Leaves the data uninitialised.
      %
      %   tmatrix = Tmatrix(data, ...)
      %   Initializes the data with the matrix `data`.
      %
      % Parameters
      %   - data (NxM numeric | cell) -- The T-matrix data.
      %     Typically a sparse or full matrix.  Data must be empty or
      %     valid T-matrix size.  If cell, describes T-matrix array
      %     and elements must be a cell array of NxM matrices.
      %
      % Optional named arguments
      %   - type (enum) -- Type of T-matrix.  Must be 'internal',
      %     'scattered' or 'total'.  Default: ``'scattered'``.
      %
      %   - index_relative (numeric) -- Relative refractive index
      %     of particle.  Default: ``1.0``.

      p = inputParser;
      p.addOptional('data', []);
      p.addParameter('type', 'scattered');
      p.addParameter('index_relative', 1.0);
      p.parse(varargin{:});

      tm = tm@ott.tmatrix.Tmatrix(p.Results.data, ...
          'type', p.Results.type);
      [tm.index_relative] = deal(p.Results.index_relative);
    end
  end

  methods % Getters/setters
    function tmatrix = set.index_relative(tmatrix, val)
      assert(isnumeric(val) && isscalar(val), ...
          'relative refractive index must be numeric scalar');
      tmatrix.index_relative = val;
    end
  end
end
