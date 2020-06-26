classdef Pointmatch < ott.Tmatrix ...
    & ott.scat.utils.RelativeMediumProperty
% Constructs a T-matrix using the point matching method
%
% The point matching method is described in
%
%   T. A. Nieminen, H. Rubinsztein-Dunlop, N. R. Heckenberg
%   JQSRT 79-80, 1019-1029 (2003), 10.1016/S0022-4073(02)00336-9
%
% Supports homogeneous isotropic materials.
% This implementation includes symmetry optimisations for rotationally
% symmetric and mirror symmetric particles.
%
% Properties
%   - relativeMedium  -- Particle relative medium
%   - rtp             -- Locations used for point matching
%   - nxyz            -- Surface normals at rtp-locations
%   - zRotSymmetry    -- Z-axis rotational symmetry (1 = no symmetry)
%   - xySymmetry      -- XY mirror symmetry (logical)
%
% Static methods
%   - DefaultProgressCallback   -- Default progress call-back method
%   - FromStarShape   -- Construct T-matrix from star shape
%
% This class is based on tmatrix_pm.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    rtp         % Location for point matching (3xN numeric)
    nxyz        % Surface normals at rtp-locations (3xN numeric)
    zRotSymmetry % Z-axis rotational symmetry (1 = no symmetry)
    xySymmetry  % XY mirror symmetry (logical)
  end

  methods (Static)
    function DefaultProgressCallback(data)
      % Default progress callback for Pointmatch
      %
      % Prints the progress to the terminal.
      %
      % Usage
      %   DefaultProgressCallback(data)
      %
      % Parameters
      %   - data (struct) -- Structure with three fields: stage
      %     (either 'setup' or 'inv'), iteration (numeric), and
      %     total (numeric).

      switch data.stage
        case 'setup'
          disp(['Setup: ' num2str(data.iteration) ...
              ' / ' num2str(data.total)]);
        case 'inv'
          disp(['Inversion: ' num2str(data.iteration) ...
              ' / ' num2str(data.total)]);
        otherwise
          error('Unknown stage');
      end
    end

    function varargout = FromStarShape(varargin)
      % Construct T-matrix using point matching from a star shaped object
      %
      % Usage
      %   tmatrix = Pointmatch(shape, relativeMedium, ...)
      %   Calculate external T-matrix.
      %
      %   [external, internal] = Pointmatch(shape, relativeMedium, ...)
      %   Calculate external and internal T-matrices.
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- A star-shaped object describing
      %     the geometry (must have a valid ``starRadii`` method).
      %
      %   - relativeMedium (ott.beam.medium.RelativeMedium) -- The relative
      %     medium describing the particle's material.
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Used to convert `shape` input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., `shape` is already in relative units).
      %
      %   - Nmax (numeric) -- Size of the VSWF expansion used for the
      %     T-matrix calculation.  In some cases it can be reduced
      %     after construction.
      %     Default: ``ott.utis.ka2nmax(2*pi*shape.maxRadius)`` (may
      %     need different values to give convergence for some shapes).
      %
      %   - angulargrid ({theta, phi}) -- Angular grid of points for
      %     calculation of radii.  Default is equally spaced angles with
      %     the number of points determined by Nmax.
      %
      % Additional parameters are passed to the class constructor.

      p = inputParser;
      p.addRequired('shape', @isnumeric);
      p.addRequired('relativeMedium', ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addParameter('wavelength', 1.0);
      p.addParameter('angulargrid', []);
      p.addParameter('Nmax', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get shape
      shape = p.Results.shape ./ p.Results.wavelength;

      % Get or calculate Nmax
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        Nmax = ott.utils.ka2nmax(2*pi*shape.maxRadius);
      else
        assert(isnumeric(Nmax) && isscalar(Nmax) && Nmax > 0, ...
            'Nmax must be positive numeric scalar');
      end

      % Get or calculate angular grid
      angulargrid = p.Results.angulargrid;
      if isempty(angulargrid)

        if shape.zRotSymmetry == 0
          ntheta = 4*(Nmax + 2);
          nphi = 1;
        else
          ntheta = 2*(Nmax + 2);
          nphi = 3*(Nmax + 2)+1;
        end

        theta = ((1:ntheta)-0.5)/ntheta * pi;
        phi = ((1:nphi)-1)/nphi * 2*pi;

        [theta, phi] = meshgrid(theta, phi);

      else
        assert(iscell(angulargrid) && numel(angulargrid) == 2, ...
            'angulargrid must be 2 element cell array');
        theta = angulargrid{1};
        phi = angulargrid{2};

        assert(isnumeric(theta), 'theta must be numeric');
        assert(isnumeric(phi), 'phi must be numeric');
        assert(numel(theta) == numel(phi), ...
            'number of theta and phi points must match');
      end

      % Calculate shape radii and normals
      radii = shape.starRadii(theta, phi);
      rtp = [radii(:), theta(:), phi(:)].';

      % Calculate surface normals
      nxyz = shape.normalsRtp(rtp);

      % Construct T-matrices
      [varargout{1:nargout}] = ott.scat.vswf.Pointmatch(...
          rtp, nxyz, p.Results.relativeMedium, 'Nmax', Nmax, ...
          'zRotSymmetry', shape.zRotSymmetry, ...
          'xySymmetry', shape.xySymmetry);
    end
  end

  methods
    function [Texternal, Tinternal] = Pointmatch(varargin)
      % Calculates T-matrix using the point matching method.
      %
      % Usage
      %   tmatrix = Pointmatch(rtp, nxyz, relativeMedium, ...)
      %   Calculate external T-matrix.
      %
      %   [external, internal] = Pointmatch(rtp, nxyz, relativeMedium, ...)
      %   Calculate external and internal T-matrices.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Locations of surface points to point-match.
      %
      %   - nxyz (3xN numeric) -- Normals at surface locations.
      %
      %   - relativeMedium (ott.beam.medium.RelativeMedium) -- The relative
      %     medium describing the particle's material.
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Used to convert `rtp` input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., `rtp(1, :)` is already in relative units).
      %
      %   - Nmax (numeric) -- Size of the VSWF expansion used for the
      %     T-matrix calculation.  In some cases it can be reduced
      %     after construction.
      %     Default: ``ott.utis.ka2nmax(2*pi*max(rtp(1, :)))`` (may
      %     need different values to give convergence for some shapes).
      %
      %   - zRotSymmetry (numeric) -- Degree of rotational symmetry
      %     about the z-axis.  Default: ``1`` (no symmetry).
      %
      %   - xySymmetry (logical) -- If the particle is mirror symmetry
      %     about the xy-plane.  Default: ``false``.
      %
      %   - progress (function_handle) -- Function to call for progress
      %     updates during method evaluation.  Takes one argument, see
      %     :meth:`DefaultProgressCallback` for more information.
      %     Default: ``[]`` (for Nmax < 20) and
      %     ``@DefaultProgressCallback`` (otherwise).

      p = inputParser;
      p.addRequired('rtp', @isnumeric);
      p.addRequired('nxyz', @isnumeric);
      p.addRequired('relativeMedium', ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addParameter('wavelength', 1.0);
      p.addParameter('Nmax', []);
      p.addParameter('progress', []);
      p.addParameter('zRotSymmetry', 1);
      p.addParameter('xySymmetry', false);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      Texternal = Texternal@ott.scat.vswf.Tmatrix(unmatched{:});
      Texternal.zRotSymmetry = p.Results.zRotSymmetry;
      Texternal.xySymmetry = p.Results.xySymmetry;
      Texternal.nxyz = p.Results.nxyz;
      Texternal.relativeMedium = p.Results.relativeMedium;
      Texternal.rtp = p.Results.rtp;
      Texternal.rtp(1, :) = Texternal.rtp(1, :) ./ p.Results.wavelength;

      % Get or calculate Nmax
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        Nmax = ott.utils.ka2nmax(2*pi*max(Texternal.rtp(1, :)));
      else
        assert(isnumeric(Nmax) && isscalar(Nmax) && Nmax > 0, ...
            'Nmax must be positive numeric scalar');
      end

      % Handle default argument for progress callback
      progress_cb = p.Results.progress;
      if isempty(progress_cb)
        if Nmax > 20
          progress_cb = @ott.scat.vswf.Pointmatch.DefaultProgressCallback;
        else
          progress_cb = @(x) [];
        end
      end

      % Calculate coefficient and incident wave matrices
      [coeff_matrix, incident_wave_matrix] = Texternal.setup(...
          Nmax, progress_callback);

      total_orders = ott.utils.combined_index(Nmax, Nmax);

      import ott.utils.*

      % Generate T-matrix

      if Texternal.zRotSymmetry == 0

        % Infinite rotational symmetry

        T = sparse(2*total_orders,2*total_orders);
        T2 = sparse(2*total_orders,2*total_orders);

        [n, m] = ott.utils.combined_index((1:total_orders).');

        for mi = -Nmax:Nmax

          % Scattered modes have the same m as incident modes
          axial_modes = m == mi;
          modes = [ axial_modes; axial_modes ];

          % This is for future, using different T and T2 size
          emodes = modes;
          omodes = modes;
          iomodes = [ omodes; emodes ];

          if Texternal.xySymmetry

            % Correct the incident modes to include even/odd modes
            even_modes = logical(mod(n + mi, 2));
            imodes_evn = modes & [ even_modes; ~even_modes ];
            imodes_odd = modes & [ ~even_modes; even_modes ];

            % Correct the outgoing modes to include even/odd modes
            even_modes = logical(mod(n + m, 2));
            omodes_evn = omodes & [ even_modes; ~even_modes ];
            omodes_odd = omodes & [ ~even_modes; even_modes ];
            emodes_evn = emodes & [ even_modes; ~even_modes ];
            emodes_odd = emodes & [ ~even_modes; even_modes ];
            iomodes_evn = [ omodes_evn; emodes_evn ];
            iomodes_odd = [ omodes_odd; emodes_odd ];

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_evn);
            Tcol = coeff_matrix(:, iomodes_evn) \ incident_wave_vectors;
            T(omodes_evn, imodes_evn) = Tcol(1:sum(omodes_evn), :);
            T2(emodes_evn, imodes_evn) = Tcol((1+sum(omodes_evn)):end, :);

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_odd);
            Tcol = coeff_matrix(:, iomodes_odd) \ incident_wave_vectors;
            T(omodes_odd, imodes_odd) = Tcol(1:sum(omodes_odd), :);
            T2(emodes_odd, imodes_odd) = Tcol((1+sum(omodes_odd)):end, :);

          else

            % Scatter the modes
            incident_wave_vectors = incident_wave_matrix(:, modes);
            Tcol = coeff_matrix(:, iomodes) \ incident_wave_vectors;
            T(omodes, modes) = Tcol(1:sum(omodes), :);
            T2(emodes, modes) = Tcol((1+sum(omodes)):end, :);

          end

          % Output progress
          progress_callback(struct('stage', 'inv', ...
              'iteration', mi, 'total', Nmax));

        end

      elseif Texternal.zRotSymmetry ~= 1

        % Discrete rotational symmetry

        T = sparse(2*total_orders,2*total_orders);
        T2 = sparse(2*total_orders,2*total_orders);

        [n, m] = ott.utils.combined_index((1:total_orders).');

        for mi = -Nmax:Nmax

          % Calculate which modes preseve symmetry, m = +/- ip
          axial_modes = mod(m - mi, Texternal.zRotSymmetry) == 0;
          incm_modes = m == mi;
          modes = [ axial_modes; axial_modes ];
          imodes = [ incm_modes; incm_modes ];

          % This is for future, using different T and T2 size
          emodes = modes;
          omodes = modes;
          iomodes = [ omodes; emodes ];

          if Texternal.xySymmetry

            % Correct the incident modes to include even/odd modes
            even_modes = logical(mod(n + mi, 2));
            imodes_evn = modes & [ even_modes; ~even_modes ];
            imodes_odd = modes & [ ~even_modes; even_modes ];

            % Correct the outgoing modes to include even/odd modes
            even_modes = logical(mod(n + m, 2));
            omodes_evn = omodes & [ even_modes; ~even_modes ];
            omodes_odd = omodes & [ ~even_modes; even_modes ];
            emodes_evn = emodes & [ even_modes; ~even_modes ];
            emodes_odd = emodes & [ ~even_modes; even_modes ];
            iomodes_evn = [ omodes_evn; emodes_evn ];
            iomodes_odd = [ omodes_odd; emodes_odd ];

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_evn);
            Tcol = coeff_matrix(:, iomodes_evn) \ incident_wave_vectors;
            T(omodes_evn, imodes_evn) = Tcol(1:sum(omodes_evn), :);
            T2(emodes_evn, imodes_evn) = Tcol((1+sum(omodes_evn)):end, :);

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_odd);
            Tcol = coeff_matrix(:, iomodes_odd) \ incident_wave_vectors;
            T(omodes_odd, imodes_odd) = Tcol(1:sum(omodes_odd), :);
            T2(emodes_odd, imodes_odd) = Tcol((1+sum(omodes_odd)):end, :);

          else

            % Scatter the modes
            incident_wave_vectors = incident_wave_matrix(:, imodes);
            Tcol = coeff_matrix(:, iomodes) \ incident_wave_vectors;
            T(omodes, imodes) = Tcol(1:sum(omodes), :);
            T2(emodes, imodes) = Tcol((1+sum(omodes)):end, :);

          end

          % Output progress
          progress_callback(struct('stage', 'inv', ...
              'iteration', mi, 'total', Nmax));

        end

      else

        if Texternal.xySymmetry

          % Only mirror symmetry
          % Parity is conserved, even modes go to even modes, etc.
          % Reference: https://arxiv.org/pdf/physics/0702045.pdf

          T = sparse(2*total_orders,2*total_orders);
          T2 = sparse(2*total_orders,2*total_orders);

          [n, m] = ott.utils.combined_index((1:total_orders).');
          even_modes = logical(mod(n + m, 2));
          modes = [ even_modes; ~even_modes ];

          % This is for future, using different T and T2 size
          imodes = modes;
          omodes = modes;
          iomodes = [ omodes; imodes ];

          % Solve for the even scattered modes
          incident_wave_vectors = incident_wave_matrix(:, modes);
          Tcol = coeff_matrix(:, iomodes) \ incident_wave_vectors;
          T(omodes, modes) = Tcol(1:sum(omodes), :);
          T2(imodes, modes) = Tcol((1+sum(omodes)):end, :);

          % Solve for the odd scattered modes
          incident_wave_vectors = incident_wave_matrix(:, ~modes);
          Tcol = coeff_matrix(:, ~iomodes) \ incident_wave_vectors;
          T(~omodes, ~modes) = Tcol(1:sum(~omodes),:);
          T2(~imodes, ~modes) = Tcol((1+sum(~omodes)):end,:);

        else

          % No rotational or mirror symmetry

          Tcol = coeff_matrix \ incident_wave_matrix;
          T = Tcol(1:2*total_orders,:);
          T2 = Tcol((1+2*total_orders):4*total_orders,:);

        end
      end

      Texternal.data = T;
      Texternal = Texternal.setType('scattered');

      if nargout == 2
        Tinternal = Texternal;
        Tinternal.data = T2;
        Tinternal = Tinternal.setType('internal');
      end
    end
  end

  methods (Access=protected)

    function [coeff_matrix, incident_wave_matrix] = setup(tmatrix, ...
        Nmax, progress_callback)
      % Calculate the coefficient and incident wave matrices

      normals = tmatrix.nxyz;
      npoints = size(tmatrix.rtp, 2);
      total_orders = ott.utils.combined_index(Nmax, Nmax);

      r = tmatrix.rtp(1, :).';
      theta = tmatrix.rtp(2, :).';
      phi = tmatrix.rtp(3, :).';

      % 3 vector components at each point, c/d,p/q coefficient per order
      coeff_matrix = zeros(6*npoints,4*total_orders);
      incident_wave_matrix = zeros(6*npoints,2*total_orders);

      k_relative = tmatrix.relativeMedium.wavenumber;

      import ott.utils.vswf;
      import ott.utils.perpcomponent;

      for n = 1:Nmax
        for m = -n:n

          % INCIDENT-SCATTERED
          [M1,N1,~,~,M2,N2] = vswf(n,m,tmatrix.k_medium*r,theta,phi);
          [M3,N3] = vswf(n,m,tmatrix.k_particle*r,theta,phi,3);

          ci = ott.utils.combined_index(n,m);

          M1 = perpcomponent(M1,normals);
          N1 = perpcomponent(N1,normals);
          M2 = perpcomponent(M2,normals);
          N2 = perpcomponent(N2,normals);
          M3 = perpcomponent(M3,normals);
          N3 = perpcomponent(N3,normals);
          M1 = M1(:);
          N1 = N1(:);
          M2 = M2(:);
          N2 = N2(:);
          M3 = M3(:);
          N3 = N3(:);

          % 1 is outgoing field, 3 is particle field, 2 is incoming field
          coeff_matrix(:,ci) = - [ M1; N1 ];
          coeff_matrix(:,ci+total_orders) = - [ N1; M1 ];
          coeff_matrix(:,ci+2*total_orders) = [ M3; k_relative*N3 ];
          coeff_matrix(:,ci+3*total_orders) = [ N3; k_relative*M3 ];

          incident_wave_matrix(:,ci) = [ M2; N2 ];
          incident_wave_matrix(:,ci+total_orders) = [ N2; M2 ];

        end

        % Output progress
        progress_callback(struct('stage', 'setup', ...
            'iteration', ci, 'total', total_orders));
      end
    end
  end

  methods % Getters/setter
    function tmatrix = set.rtp(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'rtp must be a 3xN matrix');
      tmatrix.rtp = val;
    end

    function tmatrix = set.nxyz(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'nxyz must be a 3xN matrix');
      tmatrix.nxyz = val;
    end

    function tmatrix = set.zRotSymmetry(tmatrix, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0 ...
          && round(val) == val, ...
          'zRotSymmetry must be numeric scalar integer');
      tmatrix.zRotSymmetry = val;
    end

    function tmatrix = set.xySymmetry(tmatrix, val)
      assert(islogical(val) && isscalar(val), ...
          'xySymmetry must be scalar logical');
      tmatrix.xySymmetry = val;
    end
  end
end
