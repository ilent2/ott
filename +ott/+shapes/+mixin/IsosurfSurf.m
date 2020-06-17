classdef IsosurfSurf
% Declares surface methods and casts that use an isosurface
%
% Methods
%   - surf
%   - surfPoints
%   - normalXyzInternal
%   - normalRtpInternal
%
% Supported casts
%   - PatchMesh   -- Uses isosurface
%
% All these methods raise an error when used.  These shapes do
% not support surface specific methods (for some shapes, this may
% change in future).

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function varargout = surfPoints(shape, varargin)
      % Cast to PatchMesh and call surfPoints
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)

      shape = ott.shapes.PatchMesh(shape);
      [varargout{1:nargout}] = shape.surfPoints(varargin{:});
    end

    function shape = ott.shapes.PatchMesh(shape, varargin)
      % Convert shape to a PatchMesh using default isosurface

      FV = shape.isosurface('origin', 'local');
      shape = ott.shapes.PatchMesh(FV.vertices.', FV.faces.', ...
          'position', shape.position, 'rotation', shape.rotation);
    end
  end

  methods (Hidden)
    function varargout = normalsXyzInternal(shape, varargin)
      shape = ott.shapes.PatchMesh(shape);
      [varargout{1:nargout}] = shape.normalsXyzInternal(varargin{:});
    end

    function varargout = normalsRtpInternal(shape, varargin)
      shape = ott.shapes.PatchMesh(shape);
      [varargout{1:nargout}] = shape.normalsRtpInternal(varargin{:});
    end

    function S = surfInternal(shape)
      % Cast to PatchMesh and call surf.
      %
      % Adds an {'EdgeColor', 'none'} field to the structure.
      %
      % Usage
      %   p = shape.surfInternal(...)

      shape = ott.shapes.PatchMesh(shape);
      S = shape.surfInternal();
      S.EdgeColor = 'none';
    end
  end
end
