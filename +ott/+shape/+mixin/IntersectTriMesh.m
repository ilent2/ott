classdef IntersectTriMesh
% Implement intersect by casting to a TriangularMesh
%
% Methods
%   - intersectAllInternal
%   - intersectInternal

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Hidden)
    function varargout = intersectAllInternal(shape, varargin)
      % Cast to TriangularMesh and calculate intersection

      shape = ott.shape.TriangularMesh(shape);
      [varargout{1:nargout}] = shape.intersectAllInternal(varargin{:});
    end

    function varargout = intersectInternal(shape, varargin)
      % Cast to TriangularMesh and calculate intersection

      shape = ott.shape.TriangularMesh(shape);
      [varargout{1:nargout}] = shape.intersectInternal(varargin{:});
    end
  end
end
