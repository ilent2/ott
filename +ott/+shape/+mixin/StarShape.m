classdef StarShape < ott.shape.mixin.CoordsSph
% Define methods for star shaped particles
% Inherits from :class:`CoordsSph`
%
% Properties
%   - starShaped -- Constant (true)
%
% Methods
%   - surfInternal      -- Generate a visualisation by casting to a PatchMesh
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
    % starRadii   % Not actually abstract, implemented to conflict with Shape
  end

  methods
    function shape = ott.shape.PatchMesh(shape, varargin)
      % Cast shape to a PatchMesh
      %
      % Usage
      %   shape = ott.shape.PatchMesh(shape, ...)
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

      % Calculate radii
      R = shape.starRadii(T(1:end-1, :), P(1:end-1, :));

      % Ensure ends and edge locations match
      % Reduces work and makes FromSurfMatrix behave nicer
      R(:, 1) = R(1, 1);
      R(:, end) = R(1, end);
      T(end, :) = T(1, :);
      P(end, :) = P(1, :);
      R(end+1, :) = R(1, :);

      % Convert to Cartesian coordinates
      [X, Y, Z] = ott.utils.rtp2xyz(R, T, P);

      shape = ott.shape.PatchMesh.FromSurfMatrix(X, Y, Z, ...
          'position', shape.position, 'rotation', shape.rotation);
    end

    function R = starRadii(shape, theta, phi)
      % Method implemented to cause conflict with method from Shape
      % This method should be implemented in the sub-class
      error('Method should be overloaded in sub-class');
    end
  end

  methods (Hidden)
    function S = surfInternal(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Converts the shape to a PatchMesh and calls surf.
      %
      % Usage
      %   S = shape.surfInternal(...)
      %
      % Optional named parameters
      %   - resolution ([ntheta, nphi]) -- Number of faces in theta/phi
      %     directions. Default: ``[20, 20]``.
      %
      % Additional named parameters are passed to PatchMesh.surfInternal.

      p = inputParser;
      p.addParameter('resolution', [20, 20]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      shape = ott.shape.PatchMesh(shape, 'resolution', p.Results.resolution);
      S = shape.surfInternal(unmatched{:});
    end

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
