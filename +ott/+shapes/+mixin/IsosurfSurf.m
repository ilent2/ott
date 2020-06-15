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
    function varargout = surf(shape, varargin)
      % Cast to PatchMesh and call surf
      %
      % Usage
      %   p = shape.surf(...)
      %
      % Optional parameters
      %   - surfOptions (cell) -- Optional parameters to pass to patch.
      %     Default: ``{'EdgeColor', 'none'}
      %
      % All other parameters are passed to PatchMesh.

      p = inputParser;
      p.addParameter('surfOptions', {'EdgeColor', 'none'});
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shapes.PatchMesh(shape);
      [varargout{1:nargout}] = shape.surf(...
          'surfOptions', p.Results.surfOptions, unmatched{:});

      % Change the lighting
      camlight; lighting phong;
    end

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
  end
end
