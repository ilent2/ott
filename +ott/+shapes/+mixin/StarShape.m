classdef StarShape < ott.shapes.mixin.CoordsSph
% Define methods for star shaped particles
% Inherits from :class:`CoordsSph`
%
% Properties
%   - starShaped -- Constant (true)
%
% Methods
%   - surf      -- Generate a visualisation by casting to a PatchMesh
%   - insideRtpInternal -- Uses starRadii to determine if inside
%
% Abstract methods
%   - starRadii -- Calculate the radii of the star shaped particle
%   - normalsRtpInternal -- Surface normals
%
% Supported casts
%   - PatchMesh

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    starShaped
  end

  methods (Abstract)
    starRadii
  end

  methods
    function surf(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Converts the shape to a PatchMesh and calls surf.
      %
      % Usage
      %   shape.surf(...)
      %
      % Optional named parameters
      %   - resolution ([ntheta, nphi]) -- Number of faces in theta/phi
      %     directions. Default: ``[20, 20]``.
      %
      % Additional named parameters are passed to PatchMesh.surf.

      p = inputParser;
      p.addParameter('resolution', [20, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shapes.PatchMesh(shape, 'resolution', p.Results.resolution);
      shape.surf(unmatched{:});
    end

    function shape = ott.shapes.PatchMesh(shape, varargin)
      % Cast shape to a PatchMesh
      %
      % Usage
      %   shape = ott.shapes.PatchMesh(shape, ...)
      %
      % Optional named parameters
      %   - resolution ([ntheta, nphi]) -- Number of faces in theta/phi
      %     directions.  Default: ``[20, 20]``.
      %
      % Additional named parameters are passed to constructor.

      p = inputParser;
      p.addParameter('resolution', [20, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      res = p.Results.resolution;
      theta = linspace(0, pi, res(1)+1);
      phi = linspace(0, 2*pi, res(2)+1);
      [T, P] = meshgrid(theta, phi);

      R = shape.starRadii(T, P);
      [X, Y, Z] = ott.utils.rtp2xyz(R, T, P);

      shape = ott.shapes.PatchMesh.FromSurfMatrix(X, Y, Z, ...
          'position', shape.position, 'rotation', shape.rotation);
    end
  end

  methods (Hidden)
    function b = insideRtpInternal(shape, rtp)
      % Determine if point is inside shape
      b = rtp(1, :) <= shape.starRadii(rtp(2, :), rtp(3, :));
    end
  end

  methods % Getters/setters
    function s = get.starShaped(shape)
      s = true;
    end
  end
end
