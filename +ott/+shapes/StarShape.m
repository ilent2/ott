classdef StarShape < ott.shapes.Shape
%StarShape abstract class for star shaped particles
%
% Abstract methods:
%   radii           Calculates the particle radii for angular coordinates
%   normals         Calculates the particle normals for angular coorindates
%   axialSymmetry   Returns x, y, z rotational symmetry (0 for infinite)
%   mirrorSymmetry  Returns x, y, z mirror symmetry
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
  end

  properties (Dependent)
    maxRadius       % Maximum particle radius (useful for Nmax calculation)
  end

  methods (Abstract)
    radii(shape, theta, phi);
    normals(shape, theta, phi);
    axialSymmetry(shape);
    mirrorSymmetry(shape);

    get_maxRadius(shape, varargin);
  end

  methods
    function r = get.maxRadius(shape)
      % Get the particle max radius
      r = shape.get_maxRadius();
    end

    function varargout = locations(shape, theta, phi)
      % LOCATIONS calculates Cartessian coordinate locations for points

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      [varargout{1:nargout}] = ott.utils.rtp2xyz(...
          shape.radii(theta, phi), theta, phi);
    end

    function varargout = surf(shape, varargin)
      % SURF generate a visualisation of the shape
      %
      % SURF() displays a visualisation of the shape in the current figure.
      %
      % [X, Y, Z] = surf() calculates the coordinates and arranges them
      % in a grid.

      % Create a grid of points to plot
      sz = [100, 100];
      [theta, phi] = shape.angulargrid('full', true, 'size', sz);

      % Calculate Cartesian coordinates
      [X, Y, Z] = shape.locations(theta, phi);

      % Reshape the output to match the size
      X = reshape(X, sz);
      Y = reshape(Y, sz);
      Z = reshape(Z, sz);

      % Complete the sphere (add the missing faces)
      X(:, end+1) = X(:, 1);
      Y(:, end+1) = Y(:, 1);
      Z(:, end+1) = Z(:, 1);

      % Generate the surface
      if nargout == 0
        surf(X, Y, Z, varargin{:});
      else
        varargout = { X, Y, Z };
      end
    end

    function n = normalsXyz(shape, theta, phi)
      % NORMALSXYZ calculates Cartessian normals

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      % Convert the normals to cartesian
      n = ott.utils.rtpv2xyzv(shape.normals(theta, phi), ...
          [ zeros(size(theta)), theta, phi ]);
    end

    function varargout = angulargrid(shape, varargin)
      % ANGULARGRID calculate the angular grid and radii for the shape
      %
      % This is the default function with no symmetry optimisations.
      %
      % [theta, phi] = ANGULARGRID(Nmax) gets the default angular
      % grid for the particle.
      %
      % rtp = ANGULARGRID(Nmax) or [r, theta, phi] = ANGULARGRID(Nmax)
      % calculate the radii for the locations theta, phi.
      %
      % ANGULARGRID() uses a default Nmax of 100.
      %
      % ANGULARGRID(..., 'full', full) calculates
      % an angular grid over the full sphere.
      %
      % ANGULARGRID(..., 'size', [ntheta, nphi]) uses ntheta and
      % nphi instead of Nmax for angular grid.

      p = inputParser;
      p.addOptional('Nmax', 100);
      p.addParameter('full', false);    % Not used
      p.addParameter('size', []);    % Not used
      p.parse(varargin{:});

      if isempty(p.Results.size)
        ntheta = 2*(p.Results.Nmax + 2);
        nphi = 3*(p.Results.Nmax + 2) + 1;
      else
        ntheta = p.Results.size(1);
        nphi = p.Results.size(2);
      end

      % Calculate the angular grid
      [theta, phi] = ott.utils.angulargrid(ntheta, nphi);

      if nargout == 2
        varargout{1} = theta;
        varargout{2} = phi;
      else
        % Calculate the radii
        r = shape.radii(theta, phi);
        if nargout == 1
          varargout{1} = [ r theta phi ];
        else
          varargtou{1} = r;
          varargtou{2} = theta;
          varargtou{3} = phi;
        end
      end
    end
  end
end
