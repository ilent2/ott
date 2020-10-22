classdef Mie < ott.tmatrix.Tmatrix
% Construct T-matrix with Mie scattering coefficients.
% Inherits from :class:`ott.tmatrix.Tmatrix`.
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
%   - radius          -- Radius of sphere
%   - relative_index  -- Relative refractive index
%   - relative_permeability -- Relative permeability
%
% Static methods
%   - FromShape       -- Uses ShapeMaxRadius but raises a warning
%   - ShapeVolume     -- Construct with radius set from particle volume
%   - ShapeMaxRadius  -- Construct with radius set from particle max radius
%
% See base class for additional methods/properties.
%
% This class is based on tmatrix_mie.m and tmatrix_mie_layered.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    radius            % Radius of sphere [units of wavelength]
    relative_index    % Relative refractive index
    relative_permeability % Relative permeability
  end

  methods (Static)
    function tmatrix = FromShape(shape, varargin)
      % Construct a T-matrix from a shape object.
      %
      % Uses :meth:`ShapeMaxRadius` but raises a warning if the
      % shape isn't a sphere.
      %
      % Usage
      %   tmatrix = Mie.FromShape(shape, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- The shape input.
      %
      % All other parameters passed to constructor.

      assert(numel(shape) == 1 && isa(shape, 'ott.shape.Shape'), ...
          'shape must be a single ott.shape.Shape');

      if ~isa(shape, 'ott.shape.mixin.IsSphereAbsProp') ...
          || (isa(shape, 'ott.shape.mixin.IsSphereAbsProp') ...
          && ~shape.isSphere)
        warning('ott:scat:vswf:Mie:not_a_sphere', ...
            'Shape is not a sphere, using maxRadius property');
      end

      tmatrix = ott.tmatrix.Mie.ShapeMaxRadius(shape, varargin{:});
    end

    function tmatrix = ShapeVolume(shape, varargin)
      % Construct Mie T-matrix with radius from shape volume
      %
      % Usage
      %   tmatrix = Mie.ShapeVolume(shape, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- The shape input.
      %
      % All other parameters passed to constructor.

      assert(numel(shape) == 1 && isa(shape, 'ott.shape.Shape'), ...
          'shape must be a single ott.shape.Shape');

      radius = ((3/4/pi) * shape.volume).^(1/3);
      tmatrix = ott.tmatrix.Mie(radius, varargin{:});
    end

    function tmatrix = ShapeMaxRadius(shape, varargin)
      % Construct Mie T-matrix with radius from shape max radius
      %
      % Usage
      %   tmatrix = Mie.ShapeMaxRadius(shape, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- The shape input.
      %
      % All other parameters passed to constructor.

      assert(numel(shape) == 1 && isa(shape, 'ott.shape.Shape'), ...
          'shape must be a single ott.shape.Shape');

      tmatrix = ott.tmatrix.Mie(shape.maxRadius, varargin{:});
    end
  end

  methods
    function [tmatrix, iTmatrix] = Mie(varargin)
      % Construct a new Mie T-matrix for a sphere.
      %
      % Usage
      %   tmatirx = Mie(radius, relative_index, ...)
      %   Calculate the external T-matrix (unless `internal = true`).
      %
      %   [external, internal] = Mie(radius, relative_index, ...)
      %   Calculate both the internal and external T-matrices.
      %
      % Parameters
      %   - radius (numeric) -- Radius of the sphere (in wavelength units).
      %
      %   - relative_index (numeric) -- The relative refractive index
      %     of the particle compared to the surrounding medium.
      %
      % Optional named parameters
      %   - relative_permeability (numeric) -- Relative permeability of
      %     the particle compared to surrounding medium.  Default: ``1.0``.
      %
      %   - Nmax (numeric) -- Size of the VSWF expansion used for the
      %     T-matrix calculation.
      %     Default: ``ott.utis.ka2nmax(2*pi*radius)`` (external)
      %     or ``ott.utis.ka2nmax(2*pi*radius*relative_index)`` (internal).
      %
      %   - internal (logical) -- If true, the returned T-matrix is
      %     an internal T-matrix.  Ignored for two outputs.
      %     Default: ``false``.

      p = inputParser;
      p.addOptional('radius', @isnumeric);
      p.addOptional('relative_index', 1.0, @isnumeric);
      p.addParameter('relative_permeability', 1.0, @isnumeric);
      p.addParameter('internal', false);
      p.addParameter('Nmax', []);
      p.parse(varargin{:});

      tmatrix.relative_index = p.Results.relative_index;
      tmatrix.relative_permeability = p.Results.relative_permeability;
      tmatrix.radius = p.Results.radius;

      m = tmatrix.relative_index;
      mu = tmatrix.relative_permeability;

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
            'Nmax must be positive numeric scalar');
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
      Texternal=sparse(1:2*clength, 1:2*clength, [a(indexing);b(indexing)]);

      % Internal: Calculate internal T-matrix
      d = ( h0d.*j0 - j0d.*h0 ) ./ ( m*j1.*h0d - j1d.*h0 );
      c = ( j0d.*h0 - h0d.*j0 ) ./ ( m*j1d.*h0 - h0d.*j1 );
      Tinternal=sparse(1:2*clength, 1:2*clength, [c(indexing);d(indexing)]);

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

  methods % Getters/setters
    function tmatrix = set.radius(tmatrix, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'radius must be positive numeric scalar');
      tmatrix.radius = val;
    end

    function tmatrix = set.relative_index(tmatrix, val)
      assert(isnumeric(val) && isscalar(val), ...
          'relative refractive index must be numeric scalar');
      tmatrix.relative_index = val;
    end

    function tmatrix = set.relative_permeability(tmatrix, val)
      assert(isnumeric(val) && isscalar(val), ...
          'relative permeability must be numeric scalar');
      tmatrix.relative_permeability = val;
    end
  end
end

