classdef IsosurfSurfPoints
% Declares a surfPoints method that uses isosurface
%
% Methods
%   - surfPoints

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function varargout = surfPoints(shape, varargin)
      % Convert to a patch via isosurface
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)

      FV = shape.isosurface('origin', 'local');
      shape = ott.shape.PatchMesh(FV.vertices.', FV.faces.', ...
          'position', shape.position, 'rotation', shape.rotation);

      [varargout{1:nargout}] = shape.surfPoints(varargin{:});
    end
  end
end
