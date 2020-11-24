classdef Pointmatch < ott.tmatrix.Homogeneous
% Constructs a T-matrix using the point matching method.
% Inherits from :class:`ott.tmatrix.Homogeneous`.
%
% The point matching method is described in
%
%   T. A. Nieminen, H. Rubinsztein-Dunlop, N. R. Heckenberg
%   JQSRT 79-80, 1019-1029 (2003), 10.1016/S0022-4073(02)00336-9
%
% The method can be used to construct T-matrices for star shaped particles
% with aspect ratios close to unity. Supports homogeneous isotropic materials.
% This implementation includes symmetry optimisations for rotationally
% symmetric and mirror symmetric particles.
%
% Properties
%   - index_relative  -- Relative refractive index of particle
%   - rtp             -- Locations used for point matching
%   - nrtp            -- Surface normals at rtp-locations (spherical coords.)
%   - zRotSymmetry    -- Z-axis rotational symmetry (1 = no symmetry)
%   - xySymmetry      -- XY mirror symmetry (logical)
%
% Static methods
%   - DefaultProgressCallback   -- Default progress call-back method
%   - FromShape   -- Construct T-matrix from star shape
%
% This class is based on tmatrix_pm.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    rtp           % Location for point matching (3xN numeric)
    nrtp          % Surface normals at rtp-locations (3xN numeric)
    zRotSymmetry  % Z-axis rotational symmetry (1 = no symmetry)
    xySymmetry    % XY mirror symmetry (logical)
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
          disp(['Inversion: ' num2str(data.iteration + data.total+1) ...
              ' / ' num2str(2*data.total + 1)]);
        otherwise
          error('Unknown stage');
      end
    end

    function varargout = FromShape(shape, varargin)
      % Construct T-matrix using point matching from a star shaped object
      %
      % Usage
      %   tmatrix = Pointmatch.FromShape(shape, index_relative, ...)
      %   Calculate external T-matrix.
      %
      %   [external, internal] = Pointmatch.FromShape(...)
      %   Calculate external and internal T-matrices.
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- A star-shaped object describing
      %     the geometry (must have a valid ``starRadii`` method).
      %
      %   - index_relative (numeric) -- Relative refractive index.
      %
      % Optional named parameters
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

      assert(numel(shape) == 1 && isa(shape, 'ott.shape.Shape'), ...
          'shape must be a single ott.shape.Shape');

      p = inputParser;
      p.addOptional('index_relative', 1.0, @isnumeric);
      p.addParameter('angulargrid', []);
      p.addParameter('Nmax', []);
      p.addParameter('internal', false);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get or calculate Nmax
      Nmax = ott.tmatrix.Tmatrix.getValidateNmax(...
          p.Results.Nmax, shape.maxRadius, ...
          p.Results.index_relative, p.Results.internal || nargout ~= 1);

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
      nrtp = ott.utils.xyzv2rtpv(nxyz, ott.utils.rtp2xyz(rtp));

      % Construct T-matrices
      [varargout{1:nargout}] = ott.tmatrix.Pointmatch(...
          rtp, nrtp, p.Results.index_relative, 'Nmax', Nmax, ...
          'zRotSymmetry', shape.zRotSymmetry, ...
          'xySymmetry', shape.xySymmetry, unmatched{:});
    end
  end

  methods
    function [Texternal, Tinternal] = Pointmatch(varargin)
      % Calculates T-matrix using the point matching method.
      %
      % Usage
      %   tmatrix = Pointmatch(rtp, nrtp, index_relative, ...)
      %   Calculate external T-matrix.
      %
      %   [external, internal] = Pointmatch(rtp, nrtp, index_relative, ...)
      %   Calculate external and internal T-matrices.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Locations of surface points to point-match.
      %
      %   - nrtp (3xN numeric) -- Normals at surface locations.
      %     Spherical coordinates.
      %
      %   - index_relative (numeric) -- Relative refractive index.
      %
      % Optional named parameters
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
      %   - internal (logical) -- If true, the returned T-matrix is
      %     an internal T-matrix.  Ignored for two outputs.
      %     Default: ``false``.
      %
      %   - progress (function_handle) -- Function to call for progress
      %     updates during method evaluation.  Takes one argument, see
      %     :meth:`DefaultProgressCallback` for more information.
      %     Default: ``[]`` (for Nmax < 20) and
      %     ``@DefaultProgressCallback`` (otherwise).

      p = inputParser;
      p.addRequired('rtp', @isnumeric);
      p.addRequired('nrtp', @isnumeric);
      p.addOptional('index_relative', 1.0, @isnumeric);
      p.addParameter('Nmax', []);
      p.addParameter('progress', []);
      p.addParameter('zRotSymmetry', 1);
      p.addParameter('xySymmetry', false);
      p.addParameter('internal', false);
      p.parse(varargin{:});

      Texternal.zRotSymmetry = p.Results.zRotSymmetry;
      Texternal.xySymmetry = p.Results.xySymmetry;
      Texternal.nrtp = p.Results.nrtp;
      Texternal.index_relative = p.Results.index_relative;
      Texternal.rtp = p.Results.rtp;
      Texternal.rtp(1, :) = Texternal.rtp(1, :);

      % Get or calculate Nmax
      Nmax = ott.tmatrix.Tmatrix.getValidateNmax(...
          p.Results.Nmax, max(Texternal.rtp(1, :)), ...
          Texternal.index_relative, p.Results.internal || nargout ~= 1);

      % Handle default argument for progress callback
      progress_cb = p.Results.progress;
      if isempty(progress_cb)
        if Nmax > 20
          progress_cb = @ott.tmatrix.Pointmatch.DefaultProgressCallback;
        else
          progress_cb = @(x) [];
        end
      end

      % Calculate coefficient and incident wave matrices
      [coeff_matrix, incident_wave_matrix] = Texternal.setup(...
          Nmax, progress_cb);

      total_orders = ott.utils.combined_index(Nmax, Nmax);

      import ott.utils.*

      % Generate T-matrix
      
      % Supress Matlab warnings about sparse arrays growing, we are often
      % memory bound before this realistically becomes a problem
      %#ok<*SPRIX>

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
          progress_cb(struct('stage', 'inv', ...
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
          progress_cb(struct('stage', 'inv', ...
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
      
      if nargout == 2

        Texternal.data = T;
        Texternal = Texternal.setType('scattered');
        
        Tinternal = Texternal;
        Tinternal.data = T2;
        Tinternal = Tinternal.setType('internal');
        
      else
        if p.Results.internal
          Texternal.data = T2;
          Texternal = Texternal.setType('internal');
        else
          Texternal.data = T;
          Texternal = Texternal.setType('scattered');
        end
      end
    end
  end

  methods (Hidden)
    function [coeff_matrix, incident_wave_matrix] = setup(tmatrix, ...
        Nmax, progress_callback)
      % Calculate the coefficient and incident wave matrices

      normals = tmatrix.nrtp.';
      npoints = size(tmatrix.rtp, 2);
      total_orders = ott.utils.combined_index(Nmax, Nmax);

      r = tmatrix.rtp(1, :).';
      theta = tmatrix.rtp(2, :).';
      phi = tmatrix.rtp(3, :).';

      % 3 vector components at each point, c/d,p/q coefficient per order
      coeff_matrix = zeros(6*npoints,4*total_orders);
      incident_wave_matrix = zeros(6*npoints,2*total_orders);

      k_medium = 2*pi;
      k_particle = 2*pi*tmatrix.index_relative;

      import ott.utils.vswf;

      for n = 1:Nmax
        m = -n:n;

        % INCIDENT-SCATTERED
        [M1,N1,~,~,M2,N2] = vswf(n,m,k_medium*r,theta,phi);
        [M3,N3] = vswf(n,m,k_particle*r,theta,phi,3);

        ci = ott.utils.combined_index(n,m);

        M1 = perpcomponent(M1);
        N1 = perpcomponent(N1);
        M2 = perpcomponent(M2);
        N2 = perpcomponent(N2);
        M3 = perpcomponent(M3);
        N3 = perpcomponent(N3);

        % 1 is outgoing field, 3 is particle field, 2 is incoming field
        coeff_matrix(:,ci) = - [ M1; N1 ];
        coeff_matrix(:,ci+total_orders) = - [ N1; M1 ];
        coeff_matrix(:,ci+2*total_orders) = [ M3; k_particle/k_medium*N3 ];
        coeff_matrix(:,ci+3*total_orders) = [ N3; k_particle/k_medium*M3 ];

        incident_wave_matrix(:,ci) = [ M2; N2 ];
        incident_wave_matrix(:,ci+total_orders) = [ N2; M2 ];

        % Output progress
        progress_callback(struct('stage', 'setup', ...
            'iteration', ci(end), 'total', total_orders));
      end
      
      function V = perpcomponent(V)
        % Helper to calculate perpendicular component
        V = permute(reshape(V, npoints, [], 3), [1, 3, 2]);
        V = V - sum(normals .* V, 2) .* normals;
        V = reshape(V, 3*npoints, []);
      end
    end
  end

  methods % Getters/setter
    function tmatrix = set.rtp(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'rtp must be a 3xN matrix');
      tmatrix.rtp = val;
    end

    function tmatrix = set.nrtp(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'nrtp must be a 3xN matrix');
      tmatrix.nrtp = val ./ vecnorm(val);
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

