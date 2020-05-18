classdef BscScalar < ott.beam.vswf.Bsc
% Base class for non-Bsc-array beams such as Gaussian or Bessel.
% Inherits from :class:`Bsc`.
%
% This class overrides the ArrayType behaviour of :class:`Bsc` so
% that multiple beams can't be placed in a single Bsc object.
% Instead, arrays of beams are formed using a :class:`ott.beam.Array`.
%
% This class should rarely need to be used directly, instead it is
% used as a base class for scalar beam types.
%
% Methods
%   - catInternal -- Creates a beam array instead of a single Bsc
%
% Static methods
%   - empty       -- Create an empty ott.beam.Array

  methods (Static)
    function beam = empty(varargin)
      % Construct a new empty beam array.
      %
      % Usage
      %   beam = BscScalar.empty(...)
      %
      % All parameters are passed to ott.beam.Array.empty()

      beam = ott.beam.Array.empty(varargin{:});
    end
  end

  methods
    function bsc = BscScalar(varargin)
      % Passes all arguments to :class:`Bsc`.
      %
      % This function is not indented for regular use.  See one of the
      % derived classes for usage.

      bsc = bsc@ott.beam.vswf.Bsc(varargin{:});
    end
  end

  methods (Hidden)
    function beam = catInternal(dim, varargin)
      % Concatenate two beam objects by forming a new array.

      sz = ones(1, dim+1);
      sz(dim) = numel(varargin);
      beam = ott.beam.Array('array', sz, varargin{:});
    end
  end
end
