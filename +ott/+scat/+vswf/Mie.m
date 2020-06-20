classdef Mie < ott.scat.vswf.Tmatrix
% Construct T-matrix with Mie scattering coefficients.
%
% The Mie coefficients describe the scattering of a sphere.
% They can also be used to give a reasonable estimate of the force
% for non-spherical particles when no other suitable method is available.
%
% This class supports both homogeneous dielectric (conductive and
% non-conductive) and magnetic isotropic materials.  For other materials,
% consider :class:`Dda` or :class:`MieLayered` (for layered spheres).
%
% Properties
%   - radius          -- Radius of sphere or radii of sphere layers
%   - relativeMedium  -- Relative medium describing particle material
%
% Static methods
%   - FromShape       -- Construct Mie from shape object
%
% See base class for additional methods/properties.
%
% This class is based on tmatrix_mie.m and tmatrix_mie_layered.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    radius    % (numeric) Sphere radius (relative units)
    relativeMedium  % (ott.beam.medium.Relative) Relative particle material
  end

  methods (Static)
    function tmatrix = FromShape(shape, varargin)
      % Construct a T-matrix from a shape object.
      %
      % Uses the shape maxRadius property for the sphere radius.
      % Note: For some shapes it might be better to use the average
      % radius instead.
      %
      % Usage
      %   tmatrix = Mie.FromShape(shape, ...)
      %
      % All other parameters passed to constructor.

      assert(numel(shape) == 1 && isa(shape, 'ott.shapes.Shape'), ...
          'shape must be a single ott.shapes.Shape');

      tmatrix = ott.scat.vswf.Mie(shape.maxRadius, varargin{:});
    end
  end

  methods
    function [tmatrix, iTmatrix] = Mie(varargin)
      % Construct a new Mie T-matrix for a sphere.
      %
      % Usage
      %   tmatirx = Mie(radius, relativeMedium, ...)
      %   Calculate the external T-matrix (unless `internal = true`).
      %
      %   [external, internal] = Mie(radius, relativeMedium, ...)
      %   Calculate both the internal and external T-matrices.
      %
      % Parameters
      %   - radius (numeric) -- Radius of the sphere (in relative units
      %     unless optional parameter `wavelength` is specified).
      %
      %   - relativeMedium (ott.beam.medium.RelativeMedium) -- The relative
      %     medium describing the particle's material.  Particle must
      %     be isotropic homogeneous magnetic or dielectric.
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Used to convert radius input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., radius is already in relative units).
      %
      %   - Nmax (numeric) -- Size of the VSWF expansion used for the
      %     T-matrix calculation.  Can be reduced after construction.
      %     Default: ``ott.utis.ka2nmax(2*pi*radius)`` (external)
      %     or ``ott.utis.ka2nmax(2*pi*radius*relative_index)`` (internal).
      %
      %   - internal (logical) -- If true, the returned T-matrix is
      %     an internal T-matrix.  Ignored for two outputs.
      %     Default: ``false``.

      p = inputParser;
      p.addRequired('radius', @isnumeric);
      p.addRequired('relativeMedium', ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addParameter('wavelength', 1.0);
      p.addParameter('internal', false);
      p.addParameter('Nmax', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      tmatrix = tmatrix@ott.scat.vswf.Tmatrix(unmatched{:});
      tmatrix.relativeMedium = p.Results.relativeMedium;
      tmatrix.radius = p.Results.radius ./ p.Results.wavelength;

      m = tmatrix.relativeMedium.index_relative;
      mu = tmatrix.relativeMedium.permeability_relative;

      r0 = 2*pi * tmatrix.radius;
      r1 = 2*pi * m * tmatrix.radius;

      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        if p.Results.internal && nargout == 1
          Nmax = ott.utils.ka2nmax(r1);
        else
          Nmax = ott.utils.ka2nmax(r0);
        end
      else
        assert(isnumeric(Nmax) && isscalar(Nmax) && Nmax > 0, ...
            'Nmax must be 1 or 2 element positive numeric vector');
      end

      n = 1:Nmax;
      clength = ott.utils.combined_index(Nmax, Nmax);
      indexing = ott.utils.combined_index(1:clength);

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      j0 = (sbesselj(n,r0)).';
      j1 = (sbesselj(n,r1)).';
      h0 = (sbesselh1(n,r0)).';
      j0d = (sbesselj(n-1,r0) - n.*sbesselj(n,r0)/r0).';
      j1d = (sbesselj(n-1,r1) - n.*sbesselj(n,r1)/r1).';
      h0d = (sbesselh1(n-1,r0) - n.*sbesselh1(n,r0)/r0).';

      % External: These differ from some definitions of a,b (a = -b, b = -a)
      b = -( mu*j1d.*j0 - m*j0d.*j1 ) ./ ( mu*j1d.*h0 - m*h0d.*j1 );
      a = -( mu*j0d.*j1 - m*j1d.*j0 ) ./ ( mu*h0d.*j1 - m*j1d.*h0 );
      Texternal=sparse(1:clength, 1:clength, [a(indexing);b(indexing)]);

      % Internal: Calculate internal T-matrix
      d = ( h0d.*j0 - j0d.*h0 ) ./ ( m*j1.*h0d - j1d.*h0 );
      c = ( j0d.*h0 - h0d.*j0 ) ./ ( m*j1d.*h0 - h0d.*j1 );
      Tinternal=sparse(1:clength, 1:clength, [c(indexing);d(indexing)]);

      % Assign outputs
      if nargout == 2

        iTmatrix = tmatrix;
        iTmatrix.data = Tinternal;
        iTmatrix = iTmatrix.setType('internal');

        tmatrix.data = Texternal;
        tmatrix = tmatrix.setType('scattered');

      else

        if p.Results.internal
          tmatrix.data = Tinternal;
          tmatrix = tmatrix.setType('internal');
        else
          tmatrix.data = Texternal;
          tmatrix = tmatrix.setType('scattered');
        end

      end
    end
  end

  methods (Hidden)
    function shape = getGeometry(tmatrix, wavelength)
      % Get sphere shape with radius in specified units
      shape = ott.shapes.Sphere(tmatrix.radius*wavelength, ...
          'position', tmatrix.position*wavelength, ...
          'rotation', tmatrix.rotation);
    end
  end

  methods % Getters/setters
    function tmatrix = set.radius(tmatrix, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'radius must be positive numeric scalar');
      tmatrix.radius = val;
    end

    function tmatrix = set.relativeMedium(tmatrix, val)
      assert(isa(val, 'ott.beam.medium.Relative'), ...
          'relativeMedium must be ott.beam.medium.Relative');
      tmatrix.relativeMedium = val;
    end
  end
end
