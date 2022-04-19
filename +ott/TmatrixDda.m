classdef TmatrixDda < ott.Tmatrix
% Constructs a T-matrix using discrete dipole approximation.
% Inherits from :class:`+ott.Tmatrix`.
%
% To construct a T-matrix with DDA, either the simple interface or
% the class constructor can be used.  Using the simple interface, the
% following should produce something similar to ``TmatrixMie``::
%
%   Tmatrix = ott.TmatrixDda.simple('sphere', 0.1, 'index_relative', 1.2);
%
% The DDA method requires a lot of memory to calculate the T-matrix.
% Most small desktop computers will be unable to calculate T-matrices
% for large particles (i.e., particles larger than a couple of wavelengths
% in diameter using 20 dipoles per wavelength).
% For these particles, consider using Geometric Optics
% or Finite Difference Time Domain method.
%
% See also TmatrixDda, simple.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Static, Hidden)
    function p = input_parser(varargin)
      % Helper for input parsing

      p = inputParser;

      p.addParameter('progress_callback', ...
        @ott.TmatrixDda.DefaultProgressCallback);
      p.addParameter('Nmax', []);
      p.addParameter('wavelength0', []);
      p.addParameter('spacing', []);

      p.addParameter('index_relative', []);

      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);

      p.addParameter('k_particle', []);
      p.addParameter('wavelength_particle', []);
      p.addParameter('index_particle', []);
      p.addParameter('polarizability', 'LDR');

      p.addParameter('z_mirror_symmetry', false);
      p.addParameter('z_rotational_symmetry', 1);
      p.addParameter('low_memory', false);
      p.addParameter('use_nearfield', false);
      p.addParameter('use_iterative', false);

      p.addParameter('modes', []);

      p.addParameter('verbose', true);

      % Fields to enable compatability with Tmatrix.simple
      p.addParameter('method', []);

      p.parse(varargin{:});
    end

    function DefaultProgressCallback(data)
      % Default progress callback function
      disp(['Iteration ' num2str(data.m) ' / ' num2str(max(data.mrange))]);
    end
  end

  methods (Static)
    function tmatrix = simple(shape, varargin)
      % Construct a T-matrix using DDA for simple shapes.
      %
      % Usage
      %   simple(shape, ...) constructs a new simple T-matrix for the given
      %   :class:`+ott.+shapes.Shape` object.
      %
      %   simple(name, parameters, ...) constructs a new T-matrix for the
      %   shape described by the name and parameters.
      %   For supported shape names, see :class:`+ott.+shapes.Shape.simple`.
      %
      % Optional named arguments:
      %   - spacing (numeric)         -- Spacing between dipoles.
      %     Default: ``wavelength_particle/20``
      %
      % For other named parameters, see :meth:`TmatrixDda`.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('parameters', []);
      p.parse(varargin{:});

      % Get a shape object from the inputs
      if ischar(shape) && ~isempty(p.Results.parameters)
        shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
        varargin = varargin(2:end);
      elseif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
        error('Must input either Shape object or string and parameters');
      end

      import ott.TmatrixDda;
      import ott.Tmatrix;

      % Parse remaining parameters
      p = TmatrixDda.input_parser(varargin{:});

      % Get or estimate Nmax from the inputs
      [~, k_particle] = Tmatrix.parser_wavenumber(p, 2*pi);

      % Get spacing
      spacing = p.Results.spacing;
      if isempty(spacing)
        spacing = 2*pi./k_particle./20;
      end

      % Get the symmetry of the shape
      [~,~, z_rotational_symmetry] = shape.axialSymmetry();
      [~,~, z_mirror_symmetry] = shape.mirrorSymmetry();

      % Should we avoid 0 for voxels?
      use_even = mod(z_rotational_symmetry, 2) == 0 || z_mirror_symmetry;

      % Calculate voxel locations
      voxels = shape.voxels(spacing, 'even', use_even);

      % Calculate the T-matrix using DDA
      tmatrix = TmatrixDda(voxels, varargin{:}, ...
        'spacing', spacing, ...
        'z_rotational_symmetry', z_rotational_symmetry, ...
        'z_mirror_symmetry', z_mirror_symmetry);
    end
  end

  methods
    function tmatrix = TmatrixDda(xyz, varargin)
      % Calculates T-matrix using discrete dipole approximation.
      %
      % Usage
      %   TmatrixDda(xyz, ...) calculates the T-matrix for the particle
      %   described by voxels xyz.  xyz is a 3xN matrix of coordinates
      %   for each voxel.
      %
      % The method supports homogenous and inhomogenous particles.
      % For homogeneous particles, specify the material as a scalar,
      % 3x1 vector or 3x3 polarizability matrix.  For inhomogeneous
      % particles use a N, 3xN or 3x3N vector/matrix.
      %
      % Optional named parameters
      %   - Nmax      [r,c]   Size of the T-matrix to generate.
      %     Default: ``ott.utils.ka2nmax(max_radius*k_medium)``
      %
      %   - k_medium (numeric)          -- Wavenumber in medium
      %   - wavelength_medium (numeric) -- Wavelength in medium
      %   - index_medium (numeric)      -- Refractive index in medium.
      %     Default: ``k_medium = 2*pi``
      %
      %   - k_particle (numeric)          -- Wavenumber in particle
      %   - wavelength_particle (numeric) -- Wavelength in particle
      %   - index_particle (numeric)      -- Refractive index in particle.
      %     Default: ``k_particle = 2*pi*index_relative``
      %   - polarizability (enum|numeric) -- Polarizability or method
      %     name to use to calculate from relative refractive index.
      %     Default: 'LDR'.  Supported methods: 'LDR', 'FCD', 'CM'.
      %   - index_relative (numeric)  -- Relative refractive index.
      %     Default: 1.0
      %
      %   - wavelength0 (numeric)     -- Wavelength in vacuum.
      %     Default: 1.0
      %
      %   - spacing (numeric) -- spacing for estimating Nmax and
      %     calculating the polarizability.  Only required when
      %     polarizability is non-numeric.
      %     Default: ``[]``
      %
      %   - z_mirror_symmetry (logical) -- If z-mirror symmetry should
      %     be used.  All voxels less than 0 are ignored.
      %     Default: ``false``.
      %   - z_rotational_symmetry (numeric) -- z-rotational symmetry.
      %     Degree of rotational symmetry.
      %     Objects with no rotational symetry should set this to 1.
      %     Default: ``1``.
      %   - low_memory (logical) -- If true, the DDA implementation
      %     favours additional calculations over additional memory use
      %     allowing simualtion of larger particles.  Default: false.
      %     Only applicable with z_mirror_symmetry and
      %     z_rotational_symmetry.
      %
      %   - modes (numeric) -- Specifies the modes to include in the
      %     calculated T-matrix.  Can either be a Nx2 matrix or a N
      %     element vector with the (n, m) modes or combined index modes
      %     respectively.  Default: ``[]``.
      %
      %   - use_nearfield (logical) -- If true, uses near-field point
      %     matching.  Default: `false`.
      %
      %   - use_iterative (logical) -- If true, uses an iterative solver.
      %     Default: ``false``.
      %
      %   - verbose (logical) -- Display additional information.
      %     Doesn't affect the display of the progress callback.
      %     Default: false.

      import ott.TmatrixDda;
      import ott.Tmatrix;

      tmatrix = tmatrix@ott.Tmatrix();

      % Parse inputs
      pa = TmatrixDda.input_parser(varargin{:});

      % Check for work
      if isempty(xyz)
        error('No voxels provided, no work to do');
      end

      % Tell the user some things
      if pa.Results.verbose
        disp('Starting TmatrixDda');
        disp(['Running with ' num2str(size(xyz, 2)) ' voxels']);
      end

      % Store inputs k_medium
      k_medium = ott.Tmatrix.parser_k_medium(pa, 2.0*pi);

      alpha = [];
      if isnumeric(pa.Results.polarizability)
        alpha = pa.Results.polarizability;
      end

      % Generate mask for filtering alpha, xyz and k_particle
      if pa.Results.z_mirror_symmetry
        mask = xyz(3, :) < 0;
      else
        mask = false(1, size(xyz, 2));
      end
      if pa.Results.z_rotational_symmetry == 0

        % TODO: Work out how to do infinite rotational symmetry DDA
        warning('Using 4-fold symmetry for voxels (infinite for PM)');
        mask = mask | (xyz(2, :) < 0 | xyz(1, :) < 0);

      elseif pa.Results.z_rotational_symmetry == 1
        % Nothing to do
      elseif pa.Results.z_rotational_symmetry == 2
        mask = mask | xyz(2, :) < 0;
      elseif pa.Results.z_rotational_symmetry == 4
        mask = mask | (xyz(2, :) < 0 | xyz(1, :) < 0);
      else
        error('Only z_rotational_symmetry == 0|2|4 supported for now');
      end

      % Filter xyz and alpha by symmetries
      alpha = tmatrix.remove_by_mask(alpha, mask);
      xyz(:, mask) = [];

      % Check we can allocate sufficient memory
      uV = memory;
      if ~pa.Results.low_memory && uV.MaxPossibleArrayBytes < (size(xyz, 2)*3)^2*8
        error('OTT:TmatrixDda:too_many_dipoles', ...
          ['May have too many voxels for calculation, ', ...
          'consider reducing particle size or voxel spacing']);
      end

      rtp = ott.utils.xyz2rtp(xyz.');

      % Tell the user some things
      if pa.Results.verbose && pa.Results.z_rotational_symmetry ~= 1
        disp(['Voxels with rotational symmetry: ' num2str(size(xyz, 2))]);
      end

      % Put everything in units of wavelengths
      wavelength_medium = 2*pi./k_medium;
      xyz = xyz ./ wavelength_medium;
      rtp(:, 1) = rtp(:, 1) ./ wavelength_medium;

      % Compute or get polarizability from inputs
      if isnumeric(pa.Results.polarizability)
        % Nothing more to do, already have alpha
        assert(isempty(pa.Results.index_particle), ...
            'index_particle not supported when using polarizability');
      else
        % Calculate alpha for remaining positions

        if isempty(pa.Results.spacing)
          error('spacing is needed for polarizability calculation');
        end

        % Get and filter k_particle, convert to index
        [~, k_particle] = tmatrix.parser_wavenumber(pa, 2*pi);
        k_particle = tmatrix.remove_by_mask(k_particle, mask);
        n_relative = k_particle./k_medium;

        % Get spacing in units of medium wavelength
        spacing = pa.Results.spacing;
        spacing = spacing / wavelength_medium;

        switch pa.Results.polarizability
          case 'LDR'
            alpha = ott.utils.polarizability.LDR(...
                spacing, n_relative);
          case 'FCD'
            alpha = ott.utils.polarizability.FCD(...
                spacing, n_relative);
          case 'CM'
            alpha = ott.utils.polarizability.CM(...
                spacing, n_relative);
          otherwise
            error('Unknown polarizability method name');
        end
      end

      % Get or estimate Nmax from the inputs
      Nmax = pa.Results.Nmax;
      if isempty(Nmax)

        % Get or estimate spacing (in units of wavelength)
        spacing = pa.Results.spacing;
        spacing = spacing ./ wavelength_medium;

        if isempty(spacing)
          warning('Using spacing 1/20 for Nmax estimation');
          spacing = 1./20;
        end

        maxRadius = max(rtp(:, 1)) + spacing/2;
        Nmax = ott.utils.ka2nmax(maxRadius * 2 * pi);
      end
      assert(isscalar(Nmax), 'Only scalar Nmax supported for now');

      % Calculate range of m-orders
      if ~isempty(pa.Results.modes)
        modes = pa.Results.modes;
        if size(modes, 2) == 2
          modes = ott.utils.combined_index(modes(:, 1), modes(:, 2));
        end
      else
        modes = 1:ott.utils.combined_index(Nmax, Nmax);
      end

      % Calculate the T-matrix
      tmatrix.data = tmatrix.calc_tmatrix(Nmax, xyz.', rtp, ...
          modes, alpha, pa.Results.progress_callback, ...
          pa.Results.z_mirror_symmetry, pa.Results.z_rotational_symmetry, ...
          pa.Results.low_memory, pa.Results.use_nearfield, ...
          pa.Results.use_iterative);

      % Store the type of T-matrix
      tmatrix.type = 'scattered';

    end
  end

  methods (Hidden, Static)

    function k_particle = remove_by_mask(k_particle, mask)
      % Warning: Function usage/definition may change without notice

      % Check that we have something to filter
      if isempty(k_particle)
        return;
      end

      if isscalar(k_particle)
        % Nothing to do
      elseif numel(k_particle) == 3
        % Nothing to do
      elseif numel(k_particle) == 9
        % Nothing to do
      elseif numel(k_particle) == numel(mask)
        k_particle(mask) = [];
      elseif numel(k_particle) == 3*numel(mask)
        k_particle(:, mask) = [];
      elseif numel(k_particle) == 9*numel(mask)
        mask = repelem(mask, 3);
        k_particle(:, mask) = [];
      else
        error('Unknown filter option');
      end
    end

    function M = cart2sph_mat(theta, phi)
      % Warning: Function usage/definition may change without notice

      M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
           cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
           -sin(phi) cos(phi) 0];

    end

    function M = sph2cart_mat(theta, phi)
      % Warning: Function usage/definition may change without notice

      M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
           cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
           -sin(phi) cos(phi) 0].';
    end

    function D = dipole_farfield(dipole_xyz, targets_xyz, M_dipole, k)
      % Calculate asymptotic far-field limit of dipole contributions
      %
      % Warning: Function usage/definition may change without notice

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(targets_xyz, 1) == 3, 'targets_xyz must be 3xN array');

      n_vec = targets_xyz ./ vecnorm(targets_xyz);
      k_vec = k .* n_vec;
      k_rd = sum(k_vec .* dipole_xyz, 1);

      n_targets = size(targets_xyz, 2);

      nn = reshape(n_vec, 3, 1, n_targets) .* reshape(n_vec, 1, 3, n_targets);

      % Reshape arrays for multiplication
      nn = reshape(nn, [3, 3, n_targets]);
      k_rd = reshape(k_rd, [1, 1, n_targets]);

      % Far-field limit of dipole approximation (only 1/r term)
      %
      % Based on equation 5 from
      % https://www.sciencedirect.com/science/article/pii/S0022407315002307
      F = -k^2 .* exp(-1i .* k_rd) .* (nn - eye(3));

      % Apply coordinate transformations for dipoles/targets (D = F * M)
      D = sum(reshape(F, [3, 3, 1, n_targets]) ...
          .* reshape(M_dipole, [1, size(M_dipole)]), 2);
      D = reshape(D, 3, 3, n_targets);

      % Change back to 3Nx3 result
      D = permute(D, [1, 3, 2]);
      D = reshape(D, [3*n_targets, 3]);
    end

    function D = dipole_nearfield(dipole_xyz, targets_xyz, M_dipole, k)
      % Calculate contribution from dipole to each near-field target
      %
      % dipole_xyz and targets_xyz must both be 3xN arrays
      % M_dipole and M_targets should be 3x3 and 3x3xN arrays
      % Vectorised implementation to hopefully run a little faster
      %
      % Warning: Function usage/definition may change without notice

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(targets_xyz, 1) == 3, 'targets_xyz must be 3xN array');

      n_targets = size(targets_xyz, 2);

      % Calculate distance and unit vector from dipole to targets
      r_vec = targets_xyz - dipole_xyz;
      r_jk = vecnorm(r_vec);
      r_hat = r_vec./r_jk;

      % Compute outer product of r_hat * r_hat
      % TODO: Half of these elements are redundent... Optimise?
      rr = reshape(r_hat, 3, 1, n_targets) .* reshape(r_hat, 1, 3, n_targets);

      % Reshape arrays for multiplication
      rr = reshape(rr, [3, 3, n_targets]);
      r_jk = reshape(r_jk, [1, 1, n_targets]);

      % Calculate Greens function between dipole and target location
      % This matches Eq 17 from https://doi.org/10.1364/JOSAA.18.001944
      F = -exp(1i*k*r_jk)./r_jk .* ...
        (k^2*(rr - eye(3)) + (1i*k*r_jk - 1)./r_jk.^2.*(3*rr - eye(3)));

      % Apply coordinate transformations for dipoles/targets (D = F * M)
      D = sum(reshape(F, [3, 3, 1, n_targets]) ...
          .* reshape(M_dipole, [1, size(M_dipole)]), 2);
      D = reshape(D, 3, 3, n_targets);

      % Change back to 3Nx3 result
      D = permute(D, [1, 3, 2]);
      D = reshape(D, [3*n_targets, 3]);
    end

    function F = nearfield_matrix_lowmem(theta, phi, r_near, k_medium, ...
        xyz, rtp, z_mirror, z_rotation, m, field_func)

      % Determine if we need extra memory for mirror symmetry
      if z_mirror
        z_mirror = 2;   % Yes
      else
        z_mirror = 1;   % No
      end

      % Calculate phase factor for rotational symmetry
      midx = 1:z_rotation;
      phase_factor = exp(1i*m*2*pi*(midx-1)/z_rotation);

      n_dipoles = size(rtp, 1);
      n_nfpts = length(theta);

      % Get coordinates for near-field points
      targets_xyz = ott.utils.rtp2xyz(r_near, theta, phi).';

      % Allocate space for results
      F = zeros(3*n_nfpts, 3, z_mirror, n_dipoles);

      for ii = 1:n_dipoles        % Loop over each dipole

        % Calculate cartesian to spherical conversion for each dipole
        M_cart2sph = ott.TmatrixDda.cart2sph_mat(...
            rtp(ii, 2), rtp(ii, 3));

        for jj = 1:z_rotation     % duplicates for z mirror symmetry

          % Calculate spherical to cartesian for each mirror version
          M_sph2cart = ott.TmatrixDda.sph2cart_mat(...
              rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation);

          dipole_xyz = ott.utils.rtp2xyz(rtp(ii, 1), ...
            rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation).';

          for kk = 1:z_mirror     % duplicates for z rotational symmetry

            % Apply mirror symmetry rotation
            if kk == 1
              M_dipole = M_sph2cart * M_cart2sph;
            else
              M_dipole = diag([1, 1, -1]) * M_sph2cart * M_cart2sph;
              dipole_xyz(3) = -dipole_xyz(3);
            end

            % Calculate columns of F
            F(:, :, kk, ii) = F(:, :, kk, ii) + ...
                field_func(...
                dipole_xyz, targets_xyz, M_dipole, k_medium) .* phase_factor(jj);
          end
        end
      end

      % Convert from 3xN*3*M*L*O to 3xN*3xM*L*O
      F = permute(F, [1, 2, 4, 3]);
      F = reshape(F, [3*n_nfpts, 3*n_dipoles, 1, z_mirror]);

    end

    function F = nearfield_matrix_total(theta, phi, r_near, k_medium, ...
        xyz, rtp, z_mirror, z_rotation, field_func)

      % Determine if we need extra memory for mirror symmetry
      if z_mirror
        z_mirror = 2;   % Yes
      else
        z_mirror = 1;   % No
      end

      n_dipoles = size(rtp, 1);
      n_nfpts = length(theta);

      % Get coordinates for near-field points
      targets_xyz = ott.utils.rtp2xyz(r_near, theta, phi).';

      % Allocate space for results
      F = zeros(3*n_nfpts, 3, z_mirror, z_rotation, n_dipoles);

      for ii = 1:n_dipoles        % Loop over each dipole

        % Calculate cartesian to spherical conversion for each dipole
        M_cart2sph = ott.TmatrixDda.cart2sph_mat(...
            rtp(ii, 2), rtp(ii, 3));

        for jj = 1:z_rotation     % duplicates for z mirror symmetry

          % Calculate spherical to cartesian for each mirror version
          M_sph2cart = ott.TmatrixDda.sph2cart_mat(...
              rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation);

          dipole_xyz = ott.utils.rtp2xyz(rtp(ii, 1), ...
            rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation).';

          for kk = 1:z_mirror     % duplicates for z rotational symmetry

            % Apply mirror symmetry rotation
            if kk == 1
              M_dipole = M_sph2cart * M_cart2sph;
            else
              M_dipole = diag([1, 1, -1]) * M_sph2cart * M_cart2sph;
              dipole_xyz(3) = -dipole_xyz(3);
            end

            % Calculate columns of F
            F(:, :, kk, jj, ii) = field_func(...
                dipole_xyz, targets_xyz, M_dipole, k_medium);
          end
        end
      end

      % Convert from 3xN*3*M*L*O to 3xN*3xM*L*O
      F = permute(F, [1, 2, 5, 4, 3]);
      F = reshape(F, [3*n_nfpts, 3*n_dipoles, z_rotation, z_mirror]);

    end

    function A = interaction_A_lowmem(k_medium, xyz, inv_alpha, z_mirror, z_rotation, m)

      % Determine if we need extra memory for mirror symmetry
      if z_mirror
        z_mirror = 2;   % Yes
      else
        z_mirror = 1;   % No
      end

      % Calculate phase factor for rotational symmetry
      midx = 1:z_rotation;
      phase_factor = exp(1i*m*2*pi*(midx-1)/z_rotation);

      n_dipoles = size(xyz, 1);

      rtp = ott.utils.xyz2rtp(xyz);

      % Allocate space for results
      A = zeros(3*n_dipoles, 3, z_mirror, n_dipoles);

      for ii = 1:n_dipoles        % Loop over each dipole

        % Calculate cartesian to spherical conversion for each dipole
        M_cart2sph = ott.TmatrixDda.cart2sph_mat(...
            rtp(ii, 2), rtp(ii, 3));

        for jj = 1:z_rotation     % duplicates for z mirror symmetry

          % Calculate spherical to cartesian for each mirror version
          M_sph2cart = ott.TmatrixDda.sph2cart_mat(...
              rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation);

          dipole_xyz = ott.utils.rtp2xyz(rtp(ii, 1), ...
            rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation).';

          for kk = 1:z_mirror     % duplicates for z rotational symmetry

            % Apply mirror symmetry rotation
            if kk == 1
              M_dipole = M_sph2cart * M_cart2sph;
            else
              M_dipole = diag([1, 1, -1]) * M_sph2cart * M_cart2sph;
              dipole_xyz(3) = -dipole_xyz(3);
            end

            % Calculate columns of A
            A(:, :, kk, ii) = A(:, :, kk, ii) - ...
                ott.TmatrixDda.dipole_nearfield(...
                dipole_xyz, xyz.', M_dipole, k_medium) .* phase_factor(jj);

            % Remove self-interaction terms
            if kk == 1 && jj == 1
              A((1:3) + 3*(ii-1), :, kk, ii) = inv_alpha(:, (1:3) + 3*(ii-1));
            end
          end
        end
      end

      % Convert from 3xN*3*M*L*O to 3xN*3xM*L*O
      A = permute(A, [1, 2, 4, 3]);
      A = reshape(A, [3*n_dipoles, 3*n_dipoles, 1, z_mirror]);

      % Put the inverse polarisability on the diagonal
      % Moved to main loop: R2018a this allocated twice as much memory
      %Ac = mat2cell(inv_alpha, 3, repmat(3, 1, n_dipoles));
      %A(:, :, 1, 1) = A(:, :, 1, 1) + blkdiag(Ac{:});
    end

    function A = interaction_A_total(k_medium, xyz, inv_alpha, z_mirror, z_rotation)

      % Determine if we need extra memory for mirror symmetry
      if z_mirror
        z_mirror = 2;   % Yes
      else
        z_mirror = 1;   % No
      end

      n_dipoles = size(xyz, 1);

      rtp = ott.utils.xyz2rtp(xyz);

      % Allocate space for results
      A = zeros(3*n_dipoles, 3, z_mirror, z_rotation, n_dipoles);

      for ii = 1:n_dipoles        % Loop over each dipole

        % Calculate cartesian to spherical conversion for each dipole
        M_cart2sph = ott.TmatrixDda.cart2sph_mat(...
            rtp(ii, 2), rtp(ii, 3));

        for jj = 1:z_rotation     % duplicates for z mirror symmetry

          % Calculate spherical to cartesian for each mirror version
          M_sph2cart = ott.TmatrixDda.sph2cart_mat(...
              rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation);

          dipole_xyz = ott.utils.rtp2xyz(rtp(ii, 1), ...
            rtp(ii, 2), rtp(ii, 3) + 2*pi*(jj-1)/z_rotation).';

          for kk = 1:z_mirror     % duplicates for z rotational symmetry

            % Apply mirror symmetry rotation
            if kk == 1
              M_dipole = M_sph2cart * M_cart2sph;
            else
              M_dipole = diag([1, 1, -1]) * M_sph2cart * M_cart2sph;
              dipole_xyz(3) = -dipole_xyz(3);
            end

            % Calculate columns of A
            A(:, :, kk, jj, ii) = -ott.TmatrixDda.dipole_nearfield(...
                dipole_xyz, xyz.', M_dipole, k_medium);

            % Remove self-interaction terms
            if kk == 1 && jj == 1
              A((1:3) + 3*(ii-1), :, kk, jj, ii) = zeros(3, 3);
            end
          end
        end
      end

      % Convert from 3xN*3*M*L*O to 3xN*3xM*L*O
      A = permute(A, [1, 2, 5, 4, 3]);
      A = reshape(A, [3*n_dipoles, 3*n_dipoles, z_rotation, z_mirror]);

      % Put the inverse polarisability on the diagonal
      Ac = mat2cell(inv_alpha, 3, repmat(3, 1, n_dipoles));
      A(:, :, 1, 1) = A(:, :, 1, 1) + blkdiag(Ac{:});
    end

    function [even, odd] = combine_rotsym_matrix(F, m, z_rotation)
      % Warning: Function usage/definition may change without notice

      % Combine z rotational symmetry slices
      if size(F, 3) ~= 1
        midx = 1:z_rotation;
        phase_factor = exp(1i*m*2*pi*(midx-1)/z_rotation);
        F = sum(F .* reshape(phase_factor, [1, 1, z_rotation, 1]), 3);
      end

      % Calculate even and odd parity reflections
      even = sum(F, 4);
      if size(F, 4) == 1
        odd = even;
      else
        odd = F(:, :, :, 1) - F(:, :, :, 2);
      end

    end

    function MN = calculate_farfield_modes(theta, phi, Nmax)
      % Pre-calculates Ms and Ns for all pts and all modes
      % Warning: Function usage/definition may change without notice

      total_orders = ott.utils.combined_index(Nmax,Nmax);
      npts = length(theta);

      roworder = [2:3:3*npts, 3:3:3*npts].';

      MN = zeros(3*npts,2*total_orders);
      for n = 1:Nmax
        m = -n:n;
        ci = ott.utils.combined_index(n,m);

        % Calculate fields
        [~,dtY,dpY]= ott.utils.spharm(n,m,theta,phi);
        Mnm = [dpY;-dtY] * (-1i)^(n+1)/sqrt(n*(n+1));
        Nnm = [dtY;dpY] * (-1i)^(n)/sqrt(n*(n+1));

        % Package output in [theta1; phi1; theta2; phi2; ...]
        MN(roworder,ci) = Mnm;
        MN(roworder,ci+total_orders) = Nnm;
      end

      % Convert spherical to Cartesian coordinates
      % This DDA implementation works in Cartesian
      % TODO: Perhaps its better to have a DDA/PM which works
      %   with Spherical coordinates instead? (i.e. put this
      %   in dipole_farfield instead)
      for ii = 1:numel(theta)
        MN((1:3) + (ii-1)*3, :) = ...
            ott.TmatrixDda.sph2cart_mat(theta(ii), phi(ii)) ...
            * MN((1:3) + (ii-1)*3, :);
      end
    end

    function MN = calculate_nearfield_modes(kr, theta, phi, Nmax)
      % precalculate M1 * N1's for all pts and all modes
      % Warning: Function usage/definition may change without notice

      total_orders = ott.utils.combined_index(Nmax,Nmax);
      npts = length(theta);

      MN = zeros(3*npts,2*total_orders);
      for n = 1:Nmax
        for m = -n:n
          ci = ott.utils.combined_index(n,m);
          [M1,N1] = ott.utils.vswfcart(n, m, kr, theta, phi, 'outgoing');

          % Package output in [theta1; phi2; theta2; phi2; ...]
          M1nm = M1.';
          N1nm = N1.';
          MN(:,ci) = M1nm(:);
          MN(:,ci+total_orders) = N1nm(:);
        end
      end
    end

    function inv_alpha = alpha_to_full_inv_alpha(alpha, n_dipoles)
      % Warning: Function usage/definition may change without notice

      sz = size(alpha);
      if all(sz == [1, 1])
        inv_alpha = repmat([1./alpha, 0, 0; ...
            0, 1./alpha, 0; 0, 0, 1./alpha], 1, n_dipoles);
      elseif all(sz == [3, 3])
        inv_alpha = repmat(inv(alpha), 1, n_dipoles);
      elseif all(sz == [3, 1])
        inv_alpha = repmat(diag(1.0./alpha), 1, n_dipoles);
      elseif numel(alpha) == n_dipoles
        inv_alpha = zeros(3, 3*n_dipoles);
        inv_alpha(1, 1:3:end) = 1.0./alpha;
        inv_alpha(2, 2:3:end) = 1.0./alpha;
        inv_alpha(3, 3:3:end) = 1.0./alpha;
      elseif all(sz == [3, n_dipoles])
        inv_alpha = zeros(3, 3*n_dipoles);
        inv_alpha(1, 1:3:end) = 1.0./alpha(1:3:end);
        inv_alpha(2, 2:3:end) = 1.0./alpha(2:3:end);
        inv_alpha(3, 3:3:end) = 1.0./alpha(3:3:end);
      elseif all(sz == [3, 3*n_dipoles])
        inv_alpha = zeros(3, 3*n_dipoles);
        for ii = 1:n_dipoles
          inv_alpha(:, (1:3) + 3*(ii-1)) = inv(alpha(:, (1:3) + 3*(ii-1)));
        end
      else
        error('Unsupported size of alpha');
      end
    end

    function [theta, phi] = angulargrid(Nmax, z_mirror, z_rotation)
      % Generate the angular grid for DDA point matching

      % We need at least one point for every beam shape coefficient
      % With no mirror or rotational symmetry, use the same grid as PmGauss
      ntheta = (Nmax + 1);
      nphi = 2*(Nmax + 1);

      % We can reduce the number of point around the z-axis when we
      % have z rotational symmetry since modes will only scatter to
      % other modes with a multiple of the rotational symmetry factor.
      if z_rotation > 1
        nphi = ceil(nphi ./ z_rotation);
      elseif z_rotation == 0

        % TODO: When we work out how to infinite rotational DDA
        %   we should change this value to 1, for now its 3.
        nphi = 3;
      end

      % Similarly, we can reduce the number of point in theta
      % when we have z-mirror symmetry since we match twice
      if z_mirror
        ntheta = ceil(ntheta ./ 2);
      end

      % Finally, generate the grid
      [theta, phi] = ott.utils.angulargrid(ntheta, nphi);

      % We also need to rescale our points by a similar amount
      % This avoids making the problem rank deficient
      if z_rotation > 1
        phi = phi ./ z_rotation;
      end
      if z_mirror
        theta = theta ./ 2;
      end
    end

    function data = calc_tmatrix(Nmax, xyz, rtp, modes, alpha, ...
        progress_callback, z_mirror, z_rotation, low_memory, ...
        use_nearfield, use_iterative)
      % Near-field implementation of T-matrix calculation
      %
      % This is based on DDA/T-matrix/near_field/DDA_T_NF.m
      % Warning: Function usage/definition may change without notice

      k = 2*pi;

      % Because we don't use infinite rotational symmetry everywhere,
      % we have a safe 4-fold symmetry variable.
      % TODO: Remove this when we work out infinite rotational symmetry
      z_rotation_safe = z_rotation;
      if z_rotation_safe == 0
        z_rotation_safe = 4;
      end

      % Generate grid of points over sphere
      [theta, phi] = ott.TmatrixDda.angulargrid(Nmax, z_mirror, z_rotation);

      total_orders = ott.utils.combined_index(Nmax,Nmax);

      if use_nearfield
        % Hmm, what if we have a larger sphere?
        r_near = 8; % near field radius
        assert(r_near > max(rtp(1, :)), ...
            'Particle radius too large for nearfield F');

        MN = ott.TmatrixDda.calculate_nearfield_modes(...
          k*r_near, theta, phi, Nmax);
        field_func = @ott.TmatrixDda.dipole_nearfield;
      else
        r_near = 1;   % Parameter largely ignored

        MN = ott.TmatrixDda.calculate_farfield_modes(theta, phi, Nmax)./k;
        field_func = @ott.TmatrixDda.dipole_farfield;
      end

      % Pre-calculate inv-alpha
      inv_alpha = ott.TmatrixDda.alpha_to_full_inv_alpha(alpha, size(rtp, 1));
      assert(all(isfinite(inv_alpha(:))), ...
          'singular polarizability not yet supported');

      if ~low_memory
        % When using mirror/rotational symmetry, A is still full size
        % This maybe produces a speed optimisation but uses lots of memory

        % Pre-calculate A
        A_total = ott.TmatrixDda.interaction_A_total(k, xyz, inv_alpha, ...
            z_mirror, z_rotation_safe);

        % Pre-calculate F
        F_total = ott.TmatrixDda.nearfield_matrix_total(...
            theta, phi, r_near, k, xyz, rtp, z_mirror, z_rotation_safe, ...
            field_func);

      end

      % Allocate memory for T-matrix
      data = zeros(2*total_orders, 2*total_orders);

      [nmodes, mmodes] = ott.utils.combined_index(modes(:));

      for m = unique(mmodes).'

        progress_callback(struct('m', m, 'mrange', mmodes));

        if low_memory

          % Pre-calculate A
          A_total = ott.TmatrixDda.interaction_A_lowmem(k, xyz, inv_alpha, ...
              z_mirror, z_rotation_safe, m);

          % Pre-calculate F
          F_total = ott.TmatrixDda.nearfield_matrix_lowmem(...
              theta, phi, r_near, k, xyz, rtp, z_mirror, z_rotation_safe, m, ...
              field_func);

        end

        [Feven, Fodd] = ott.TmatrixDda.combine_rotsym_matrix(F_total, ...
            m, z_rotation_safe);
        [Aeven, Aodd] = ott.TmatrixDda.combine_rotsym_matrix(A_total, ...
            m, z_rotation_safe);

        % Find which modes we should calculate
        ournmodes = nmodes(mmodes == m).';

        if ~use_iterative
          % Not supported by gmres

          % Pre-calculate all fields
          Ei_TE = zeros(3*size(rtp, 1), numel(ournmodes));
          Ei_TM = Ei_TE;
          for ii = 1:numel(ournmodes)
            [E, M] = ott.utils.vswfcart(ournmodes(ii), m, ...
                rtp(:, 1)*k, rtp(:, 2), rtp(:, 3), 'regular');
            E = E.';
            M = M.';
            Ei_TE(:, ii) = E(:);
            Ei_TM(:, ii) = M(:);
          end

          % Solve all n-modes together
          if z_mirror
            evn = mod(m+ournmodes, 2) == 0;

            E_odd = ott.TmatrixDda.solve_and_evaluate(Aodd, Fodd, ...
                [Ei_TE(:, evn), Ei_TM(:, ~evn)], use_iterative);
            E_evn = ott.TmatrixDda.solve_and_evaluate(Aeven, Feven, ...
                [Ei_TM(:, evn), Ei_TE(:, ~evn)], use_iterative);

            En_TE = zeros(size(E_odd, 1), numel(ournmodes));
            En_TM = En_TE;

            En_TE(:, evn) = E_odd(:, 1:sum(evn));
            En_TE(:, ~evn) = E_evn(:, sum(evn)+1:end);
            En_TM(:, evn) = E_evn(:, 1:sum(evn));
            En_TM(:, ~evn) = E_odd(:, sum(evn)+1:end);
          else
            % Solve both at the same time (to save time)
            E_EM = ott.TmatrixDda.solve_and_evaluate(Aeven, Feven, ...
                [Ei_TE, Ei_TM], use_iterative);
            En_TE = E_EM(:, 1:end/2);
            En_TM = E_EM(:, end/2+1:end);
          end
        end

        for n = ournmodes

          if ~use_iterative
            E_TE = En_TE(:, ournmodes == n);
            E_TM = En_TM(:, ournmodes == n);

          else

            % Calculate fields (Cartesian coordinates)
            [Ei_TE, Ei_TM] = ott.utils.vswfcart(n, m, ...
                rtp(:, 1)*k, rtp(:, 2), rtp(:, 3), 'regular');
            Ei_TE = Ei_TE.';
            Ei_TM = Ei_TM.';

            % Evaluate PM fields (Cartesian coordinates)
            if z_mirror
              if mod(m + n, 2) == 0
                E_TE = ott.TmatrixDda.solve_and_evaluate(Aodd, Fodd, Ei_TE(:), use_iterative);
                E_TM = ott.TmatrixDda.solve_and_evaluate(Aeven, Feven, Ei_TM(:), use_iterative);
              else
                E_TE = ott.TmatrixDda.solve_and_evaluate(Aeven, Feven, Ei_TE(:), use_iterative);
                E_TM = ott.TmatrixDda.solve_and_evaluate(Aodd, Fodd, Ei_TM(:), use_iterative);
              end
            else
              % Solve both at the same time (to save time)
              E_EM = ott.TmatrixDda.solve_and_evaluate(Aeven, Feven, ...
                  [Ei_TE(:), Ei_TM(:)], use_iterative);
              E_TE = E_EM(:, 1);
              E_TM = E_EM(:, 2);
            end
          end

          ci = ott.utils.combined_index(n,m);

          if z_rotation == 1
            % No rotational symmetry

            if z_mirror
              % Only mirror symmetry
              % Parity is conserved, even modes go to even modes, etc.
              % Reference: https://arxiv.org/pdf/physics/0702045.pdf

              % Find even modes
              [alln, allm] = ott.utils.combined_index((1:total_orders).');
              even_modes = logical(mod(alln + allm, 2));
              modes = [ even_modes; ~even_modes ];

              if ~logical(mod(n + m, 2))
                modes = ~modes;
              end

              pq1 = MN(:, modes) \ E_TE;
              pq2 = MN(:, ~modes) \ E_TM;

              data(modes,ci) = pq1;
              data(~modes,ci+total_orders) = pq2;

            else

              pq1 = MN\E_TE;
              pq2 = MN\E_TM;

              data(:,ci) = pq1;
              data(:,ci+total_orders) = pq2;

            end

          else

            [alln, allm] = ott.utils.combined_index((1:total_orders).');

            if z_rotation > 1
              % Calculate which modes preseve symmetry, m = +/- ip
              axial_modes = mod(allm - m, z_rotation) == 0;
            elseif z_rotation == 0
              % Modes only scatter to modes with same m
              axial_modes = allm == m;
            else
              error('Invalid z_rotation value');
            end

            modes = [axial_modes; axial_modes];

            if z_mirror

              % Correct the MN modes to include only even/odd modes
              even_modes = logical(mod(alln + allm, 2));

              if logical(mod(n + m, 2))
                modes_evn = modes & [ even_modes; ~even_modes ];
                modes_odd = modes & [ ~even_modes; even_modes ];
              else
                modes_odd = modes & [ even_modes; ~even_modes ];
                modes_evn = modes & [ ~even_modes; even_modes ];
              end

              pq1 = MN(:, modes_evn) \ E_TE;
              pq2 = MN(:, modes_odd) \ E_TM;

              data(modes_evn,ci) = pq1;
              data(modes_odd,ci+total_orders) = pq2;

            else

              pq1 = MN(:, modes) \ E_TE;
              pq2 = MN(:, modes) \ E_TM;

              data(modes,ci) = pq1;
              data(modes,ci+total_orders) = pq2;

            end
          end

        end		% for n = 1:m
      end		% for m = mrange

    end

    function E = solve_and_evaluate(A, F, E, use_iterative)
      % Solve DDA problem and evaluate fields at specified locations

      % Solve for poarizability of each dipole
      if use_iterative
        [P, ~] = gmres(A, E, 1, 1e-6, 30);
      else
        P = A\E;
      end

      % Calculate fields
      E = F * P;
    end
  end
end
