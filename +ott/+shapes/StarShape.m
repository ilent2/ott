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

    function varargout = locations(shape, theta, phi);
      % LOCATIONS calculates Cartessian coordinate locations for points

      theta = theta(:);
      phi = phi(:);
      [theta,phi] = ott.utils.matchsize(theta,phi);

      [varargout{1:nargout}] = ott.utils.rtp2xyz(...
          shape.radii(theta, phi), theta, phi);
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

    function varargout = angulargrid(shape, Nmax)
      % ANGULARGRID calculate the angular grid and radii for the shape
      %
      % This is the default function with no symmetry optimisations.
      %
      % [theta, phi] = ANGULARGRID(Nmax) gets the default angular
      % grid for the particle.
      %
      % rtp = ANGULARGRID(Nmax) or [r, theta, phi] = ANGULARGRID(Nmax)
      % calculate the radii for the locations theta, phi.

      ntheta = 2*(Nmax + 2);
      nphi = 3*(Nmax + 2) + 1;

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
