classdef MieLayered < ott.scat.vswf.Tmatrix ...
    & ott.scat.utils.RelativeMediumProperty
% Construct T-matrices for a layered sphere
%
% This class implements the first part of
%
%   "Improved recursive algorithm for light scattering by a multilayered
%   sphere", Wen Yang, Applied Optics 42(9), 2003
%
% and can be used to model layered spherical particles.
%
% Properties
%   - radii           -- Radii of sphere layers (inside to outside)
%   - relativeMedium  -- Relative medium array describing particle layers
%
% Static methods
%   - FromShape       -- Construct layered sphere from array of shapes
%
% See base class for additional methods/properties.
%
% This class is based on tmatrix_mie_layered.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    radii     % (numeric) Sphere radii (relative units)
  end

  methods (Static)
    function tmatrix = FromShape(shape, relativeMedium, varargin)
      % Construct a T-matrix from a shape array
      %
      % Uses the maxRadius property of each shape in the shape array
      % for the sphere radii.  Sorts the radii (and relative medium
      % entries) if radii aren't sorted.
      %
      % If the shapes aren't all centred, raises a warning.
      %
      % Usage
      %   tmatrix = MieLayered.FromShape(shape, relativeMedium, ...)
      %
      % All other parameters are passed to MieLayered.

      positions = [shape.position];
      if ~all(all(positions(:, 1) == positions))
          warning('ott:scat:vswf:MieLayered:not_centred', ...
              'shapes are not centred');
      end

      % Sort radii (and relativeMedium)
      radii = [shape.maxRadius];
      [radii, I] = sort(radii);
      relativeMedium = relativeMedium(I);

      tmatrix = ott.scat.vswf.MieLayered(radii, relativeMedium, varargin{:});
    end
  end

  methods
    function [tmatrix, iTmatrix] = MieLayered(varargin)
      % Construct a new Mie layered T-matrix
      %
      % Usage
      %   tmatirx = MieLayered(radii, relativeMedium, ...)
      %   Calculate the external T-matrix (unless `internal = true`).
      %
      %   [external, internal] = Mie(radii, relativeMedium, ...)
      %   Calculate both the internal and external T-matrices.
      %
      % Parameters
      %   - radii (numeric) -- Radius of sphere layers (in relative units
      %     unless optional parameter `wavelength` is specified).
      %     Radii should be in ascending order.
      %
      %   - relativeMedium (ott.beam.medium.RelativeMedium) -- The relative
      %     medium describing the particle material layers.  Particle must
      %     be isotropic homogeneous magnetic or dielectric.  The number
      %     of layers should match number of radii.
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Used to convert radius input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., radius is already in relative units).
      %
      %   - Nmax (numeric) -- Size of the VSWF expansion used for the
      %     T-matrix calculation.  Can be reduced after construction.
      %     Default: ``100`` (should be numerically stable).
      %
      %   - internal (logical) -- If true, the returned T-matrix is
      %     an internal T-matrix.  Ignored for two outputs.
      %     Default: ``false``.

      p = inputParser;
      p.addRequired('radii', @isnumeric);
      p.addRequired('relativeMedium', ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addParameter('wavelength', 1.0);
      p.addParameter('internal', false);
      p.addParameter('Nmax', 100);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      tmatrix = tmatrix@ott.scat.vswf.Tmatrix(unmatched{:});
      tmatrix.relativeMedium = p.Results.relativeMedium;
      tmatrix.radii = p.Results.radius ./ p.Results.wavelength;

      assert(numel(tmatrix.relativeMedium) == numel(tmatrix.radii), ...
          'number of radii must match number of relativeMedium');

      Nmax = p.Results.Nmax;

      k_layer = [tmatrix.relativeMedium.wavenumber, 2*pi];
      radius = [tmatrix.radii, tmatrix.radii(end)];

      k_layer=[tmatrix.wavenumber_particle,tmatrix.wavenumber_medium];
      radius=[tmatrix.radius,tmatrix.radius(end)];
      n_layer=k_layer/2/pi;

      n = 1:Nmax;

      lastElement=length(k_layer);

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      %generate all special functions first:
      if length(tmatrix.radii)>1
          [jN,jNd] = sbesselj(n,[k_layer.*radius, ...
              k_layer(2:end).*radius(1:end-1)]);
          [hN,hNd] = sbesselh1(n,[k_layer.*radius, ...
              k_layer(2:end).*radius(1:end-1)]);
      else
          [jN,jNd] = sbesselj(n,k_layer(:).*radius(:));
          [hN,hNd] = sbesselh1(n,k_layer(:).*radius(:));
      end
      jN=jN.';
      hN=hN.';
      jNd=jNd.';
      hNd=hNd.';

      d1_N=jNd./jN;
      d3_N=hNd./hN;
      r_N=jN./hN;

      ha_0=1/n_layer(1);
      hb_0=1;

      for ii=1:length(tmatrix.radii)
          ha_n=ha_0;
          hb_n=hb_0;

          m_0 = n_layer(ii);

          d1_0=d1_N(:,ii);
          d3_0=d3_N(:,ii);
          r_0=r_N(:,ii);

          if ii>1
              m_n = n_layer(ii-1);

              d1_n=d1_N(:,lastElement+ii-1);
              d3_n=d3_N(:,lastElement+ii-1);
              r_n=r_N(:,lastElement+ii-1);
          else
              m_n=1;

              d1_n=0;
              d3_n=0;
              r_n=0;
          end

          g0=m_0*ha_n-m_n*d1_n;
          g1=m_0*ha_n-m_n*d3_n;
          g0b=m_n*hb_n-m_0*d1_n;
          g1b=m_n*hb_n-m_0*d3_n;
          q_0=r_n./r_0;
          ha_0=(g1.*d1_0-q_0.*g0.*d3_0)./(g1-q_0.*g0);
          hb_0=(g1b.*d1_0-q_0.*g0b.*d3_0)./(g1b-q_0.*g0b);

      end
      m_0 = n_layer(lastElement-1);
      m_1 = n_layer(lastElement);

      m=m_0/m_1;

      d1_1=d1_N(:,lastElement);
      d3_1=d3_N(:,lastElement);
      r_1=r_N(:,lastElement);

      % %TM/TE coeffs...
      al_1=r_1.*(ha_0-m*d1_1)./(ha_0-m*d3_1);
      bl_1=r_1.*(m*hb_0-d1_1)./(m*hb_0-d3_1);

      %swap the modes to TE/TM... and negatize.
      b = -al_1;
      a = -bl_1;

      %t-matrix indices...
      clength = ott.utils.combined_index(Nmax, Nmax);
      indexing=ott.utils.combined_index(1:clength);

      Texternal = sparse(1:clength, 1:clength, [a(indexing);b(indexing)]);

      r_0=(jN(:,lastElement)./jN(:,lastElement-1));

      d = r_0.*(d3_1 - d1_1 )  ./ (ha_0-m*d3_1);
      c = r_0.*(d3_1 - d1_1 )  ./ (m*hb_0 - d3_1);

      Tinternal = sparse(1:clength, 1:clength, [c(indexing);d(indexing)]);

      if nargout == 2 || p.Results.internal
        warning('ott:scat:vswf:MieLayered::internal_coefficient', ...
            ['The internal coefficients are for the outermost layer, ' ...
             newline, ' only the real ones are only defined for each layer.']);
      end

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
      % Get spheres with radii in specified units
      shape = ott.shapes.Sphere(tmatrix.radii*wavelength, ...
          'position', tmatrix.position*wavelength, ...
          'rotation', tmatrix.rotation);
    end
  end

  methods % Getters/setters
    function tmatrix = set.radius(tmatrix, val)
      assert(isnumeric(val) && isvector(val) && issorted(val) ...
          && all(val > 0), ...
          'radii should be sorted positive numeric vector');
      tmatrix.radii = val(:).';
    end
  end
end
