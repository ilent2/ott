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
      
      if z_rotational_symmetry == 0
        warning('Using z_rotational_symmetry = 4 instead of 0 for now');
        z_rotational_symmetry = 4;
      end
      
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

      % Store inputs k_medium and k_particle
      [k_medium, k_particle] = tmatrix.parser_wavenumber(pa, 2*pi);
      
      alpha = [];
      if isnumeric(pa.Results.polarizability)
        alpha = pa.Results.polarizability;
      end
      
      % Filter xyz and k_particle for symmetries
      if pa.Results.z_mirror_symmetry
        k_particle = tmatrix.remove_by_mask(k_particle, xyz(3, :) < 0);
        alpha = tmatrix.remove_by_mask(alpha, xyz(3, :) < 0);
        xyz(:, xyz(3, :) < 0) = [];
      end
      if pa.Results.z_rotational_symmetry == 0
        error('z_rotational_symmetry == 0 not yet supported');
      elseif pa.Results.z_rotational_symmetry == 1
        % Nothing to do
      elseif pa.Results.z_rotational_symmetry == 2
        k_particle = tmatrix.remove_by_mask(k_particle, xyz(2, :) < 0);
        alpha = tmatrix.remove_by_mask(alpha, xyz(2, :) < 0);
        xyz(:, xyz(2, :) < 0) = [];
      elseif pa.Results.z_rotational_symmetry == 4
        k_particle = tmatrix.remove_by_mask(...
            k_particle, xyz(2, :) < 0 | xyz(1, :) < 0);
        alpha = tmatrix.remove_by_mask(...
            alpha, xyz(2, :) < 0 | xyz(1, :) < 0);
        xyz(:, xyz(2, :) < 0 | xyz(1, :) < 0) = [];
      else
        error('Only z_rotational_symmetry == 2|4 supported for now');
      end
      
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
        disp(['Voxels with rotationa symmetry: ' num2str(size(xyz, 2))]);
      end

      % Put everything in units of wavelengths
      wavelength_medium = 2*pi./k_medium;
      xyz = xyz ./ wavelength_medium;
      rtp(:, 1) = rtp(:, 1) ./ wavelength_medium;

      n_relative = k_particle./k_medium;

      % Compute or get polarizability from inputs
      if isnumeric(pa.Results.polarizability)
        % Nothing more to do, already have alpha
      else
        % Calculate alpha for remaining positions
        
        if isempty(pa.Results.spacing)
          error('Spacing is needed for polarizability calculation');
        end
        
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
          warning('Estimating spacing for Nmax calculation from k_particle');
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
      tmatrix.data = tmatrix.calc_nearfield(Nmax, xyz.', rtp, ...
          n_relative, modes, alpha, pa.Results.progress_callback, ...
          pa.Results.z_mirror_symmetry, pa.Results.z_rotational_symmetry, ...
          pa.Results.low_memory);

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
    
    function E = apply_cart2sph(E, cart2sph)
      % Apply coordinate transformation to E-field vector
      % Warning: Function usage/definition may change without notice
      
      for ii = 1:size(cart2sph, 2)/3
        M = cart2sph(:, (1:3) + 3*(ii-1));
        E((1:3) + 3*(ii-1)) = M * E((1:3) + 3*(ii-1));
      end
    end
    
    function D = dipole_contribution(dipole_xyz, targets_xyz, M_dipole, k)
      % Calculate contribution from dipole to each target
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
      F = exp(1i*k*r_jk)./r_jk .* ...
        (k^2*(rr - eye(3)) + (1i*k*r_jk - 1)./r_jk.^2.*(3*rr - eye(3)));
      
      % Apply coordinate transformations for dipoles/targets (D = F * M)
      D = sum(reshape(F, [3, 3, 1, n_targets]) .* reshape(M_dipole, [1, size(M_dipole)]), 2);
      D = reshape(D, 3, 3, n_targets);
      
      % Change back to 3Nx3 result
      D = permute(D, [1, 3, 2]);
      D = reshape(D, [3*n_targets, 3]);
    end
    
    function F = nearfield_matrix_lowmem(theta, phi, r_near, k_medium, ...
        xyz, rtp, z_mirror, z_rotation, m)
      
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
                ott.TmatrixDda.dipole_contribution(...
                dipole_xyz, targets_xyz, M_dipole, k_medium) .* phase_factor(jj);
          end
        end
      end
      
      % Convert from 3xN*3*M*L*O to 3xN*3xM*L*O
      F = permute(F, [1, 2, 4, 3]);
      F = reshape(F, [3*n_nfpts, 3*n_dipoles, 1, z_mirror]);
      
    end
    
    function F = nearfield_matrix_total(theta, phi, r_near, k_medium, ...
        xyz, rtp, z_mirror, z_rotation)
      
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
            F(:, :, kk, jj, ii) = ott.TmatrixDda.dipole_contribution(...
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
            A(:, :, kk, ii) = A(:, :, kk, ii) + ...
                ott.TmatrixDda.dipole_contribution(...
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
            A(:, :, kk, jj, ii) = ott.TmatrixDda.dipole_contribution(...
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

    function [F, Fodd] = nearfield_matrix(theta, phi, r_near, k, xyz, rtp, ...
          z_mirror, z_rotation)
      % Based on DDA/T-matrix/near_field/nearfield_matrix
      %
      % Warning: Function usage/definition may change without notice
      
      assert(z_rotation >= 1, ...
        'Only finite z rotationaly symmetry supported for now');

      n_dipoles = size(rtp, 1);
      n_nfpts = length(theta);
      
      xyz_E = ott.utils.rtp2xyz(r_near,theta,phi);

      F = zeros(3*n_nfpts,3*n_dipoles, z_rotation);
      Fodd = F;
      
      I = eye(3);
      
      for m = 1:z_rotation
        
        % Get voxel location for this quadrant
        if m == 1
          our_xyz = xyz;
        else
          our_xyz = ott.utils.rtp2xyz(rtp(:, 1), ...
            rtp(:, 2), rtp(:, 3) + 2*pi*(m-1)/z_rotation);
        end
        
        for n = 1:n_dipoles
          n_i = 3*(n-1)+1;   % F column index

          r_vec = xyz_E - our_xyz(n, :);
          r_jk = vecnorm(r_vec, 2, 2);
          r_hat = r_vec./r_jk;
          
          % Calculate location of mirror point
          if z_mirror
            rM_vec = xyz_E - [our_xyz(n, 1:2), -our_xyz(n, 3)];
            rM_jk = vecnorm(rM_vec, 2, 2);
            rM_hat = rM_vec./rM_jk;
          end
          
          % Calculate coordinate rotation
          if m == 1
            M = 1;
          else
            % Calculate z rotation coordinate transformation matrix
            M_cart2sph = ott.TmatrixDda.cart2sph_mat(...
                rtp(n, 2), rtp(n, 3));
            M_sph2cart = ott.TmatrixDda.sph2cart_mat(...
                rtp(n, 2), rtp(n, 3) + 2*pi*(m-1)/z_rotation);
            M = M_sph2cart * M_cart2sph;
          end
          
          for p = 1:n_nfpts
            p_i = 3*(p-1)+1;    % F row index
            
            rr = r_hat(p, :).'*r_hat(p, :);
            
            F(p_i:p_i+2, n_i:n_i+2, m) = exp(1i*k*r_jk(p, :))/r_jk(p, :)*...
              (k^2*(rr - I) + (1i*k*r_jk(p, :)-1)/r_jk(p, :)^2*(3*rr - I)) * M;
            
            % Add z-mirror-symmetry point
            if z_mirror
              Mm = diag([1, 1, -1]);
              rr = rM_hat(p, :).'*rM_hat(p, :);
              
              term = exp(1i*k*rM_jk(p, :))/rM_jk(p, :)*...
                (k^2*(rr - I) + (1i*k*rM_jk(p, :)-1)/rM_jk(p, :)^2*(3*rr - I)) * Mm * M;
              
              Fodd(p_i:p_i+2, n_i:n_i+2, m) = F(p_i:p_i+2, n_i:n_i+2, m) - term;
              F(p_i:p_i+2, n_i:n_i+2, m) = F(p_i:p_i+2, n_i:n_i+2, m) + term;
                
            end
          end
        end
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
          [M1,N1] = ott.utils.vswf(n, m, kr, theta, phi, 'outgoing');
          
          M1nm = M1.';
          N1nm = N1.';

          MN(:,ci) = -M1nm(:);
          MN(:,ci+total_orders) = -N1nm(:);
        end
      end
    end
    
    function inv_alpha = alpha_to_full_inv_alpha(alpha, n_dipoles)
      % Warning: Function usage/definition may change without notice
      
      sz = size(alpha);
      if all(sz == [1, 1])
        inv_alpha = repmat([1./alpha, 0, 0; 0, 1./alpha, 0; 0, 0, 1./alpha], 1, n_dipoles);
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

    function data = calc_nearfield(Nmax, xyz, rtp, n_rel, modes, alpha, ...
        progress_callback, z_mirror, z_rotation, low_memory)
      % Near-field implementation of T-matrix calculation
      %
      % This is based on DDA/T-matrix/near_field/DDA_T_NF.m
      % Warning: Function usage/definition may change without notice

      iterative = true;
      k = 2*pi;

      import ott.utils.combined_index;

      total_orders = combined_index(Nmax,Nmax);

      % What is this?
      npts = round(sqrt(8*total_orders));
      [theta,phi] = ott.utils.angulargrid(npts,npts);
      
      % Throw away near-field points not needed for rotational/mirror symmetry
      if z_mirror
        phi(theta > pi/2) = [];
        theta(theta > pi/2) = [];
      end
      if z_rotation == 0
        error('z_rotational_symmetry == 0 not yet supported');
      elseif z_rotation == 1
        % Nothing to do
      elseif z_rotation == 2
        theta(phi > pi) = [];
        phi(phi > pi) = [];
      elseif z_rotation == 4
        theta(phi > pi/2) = [];
        phi(phi > pi/2) = [];
      else
        error('Only z_rotational_symmetry == 2|4 supported for now');
      end

      % Hmm, what if we have a larger sphere?
      r_near = 8; % near field radius
      assert(r_near > max(rtp(1, :)), 'Particle radius too large');
      
      % Pre-calculate inv-alpha
      inv_alpha = ott.TmatrixDda.alpha_to_full_inv_alpha(alpha, size(rtp, 1));
      assert(all(isfinite(inv_alpha(:))), 'singular polarizability not yet supported');
      
      if ~low_memory
        % When using mirror/rotational symmetry, A is still full size
        % This maybe produces a speed optimisation but uses lots of memory

        % Pre-calculate A
        A_total = ott.TmatrixDda.interaction_A_total(k, rtp, inv_alpha, ...
            z_mirror, z_rotation);

        % Pre-calculate F
        F_total = ott.TmatrixDda.nearfield_matrix_total(...
            theta, phi, r_near, k, xyz, rtp, z_mirror, z_rotation);
          
      end

      MN = ott.TmatrixDda.calculate_nearfield_modes(...
        k*r_near, theta, phi, Nmax);
      
      % Pre-compute near-field Cartesian to Spherical transform
      cart2sph = zeros(3, 3*length(theta));
      for ii = 1:length(theta)
        cart2sph(:, (1:3) + (ii-1)*3) = ott.TmatrixDda.cart2sph_mat(...
            theta(ii), phi(ii));
      end

      % Allocate memory for T-matrix
      data = zeros(2*total_orders, 2*total_orders);
      
      [nmodes, mmodes] = ott.utils.combined_index(modes(:));

      for m = unique(mmodes).'

        progress_callback(struct('m', m, 'mrange', mmodes));

        if low_memory

          % Pre-calculate A
          A_total = ott.TmatrixDda.interaction_A_lowmem(k, rtp, inv_alpha, ...
              z_mirror, z_rotation, m);

          % Pre-calculate F
          F_total = ott.TmatrixDda.nearfield_matrix_lowmem(...
              theta, phi, r_near, k, xyz, rtp, z_mirror, z_rotation, m);

        end

        [Feven, Fodd] = ott.TmatrixDda.combine_rotsym_matrix(F_total, ...
            m, z_rotation);
        [Aeven, Aodd] = ott.TmatrixDda.combine_rotsym_matrix(A_total, ...
            m, z_rotation);
          
        % Find which modes we should calculate
        ournmodes = nmodes(mmodes == m).';

        for n = ournmodes
          
          % Calculate fields in Cartesian coordinates
          if z_mirror
            if mod(m + n, 2) == 0
              A_TM = Aeven;
              A_TE = Aodd;
              F_TM = Feven;
              F_TE = Fodd;
            else
              A_TE = Aeven;
              A_TM = Aodd;
              F_TE = Feven;
              F_TM = Fodd;
            end
          else
            A_TM = Aeven;
            A_TE = Aeven;
            F_TM = Feven;
            F_TE = Feven;
          end

          [Ei_TE, Ei_TM] = ott.utils.vswfcart(n, m, ...
              rtp(:, 1)*k, rtp(:, 2), rtp(:, 3), 'regular');
          Ei_TE = Ei_TE.';
          Ei_TM = Ei_TM.';

          if iterative
            [P_TE, ~] = gmres(A_TE,Ei_TE(:),1,1e-6,30);
            [P_TM, ~] = gmres(A_TM,Ei_TM(:),1,1e-6,30);
          else
            P_TE = A_TE\Ei_TE(:);
            P_TM = A_TM\Ei_TM(:);
          end

          % Add n_rel correction
          if numel(n_rel) == 1
            P_TE = P_TE * n_rel;
          elseif numel(n_rel) == 3
            P_TE = repmat(n_rel(:), [numel(P_TE)/3, 1]) .* P_TE;
          elseif numel(n_rel) == numel(P_TE)/3
            P_TE = repelem(n_rel(:), 3) .* P_TE;
          elseif numel(n_rel) == numel(P_TE)
            P_TE = n_rel(:) .* P_TE;
          elseif numel(n_rel) == 3*numel(P_TE) && size(n_rel, 1) == 3
            for ii = 1:(length(P_TE)/3)
              P_TE((1:3) + (ii-1)*3) = ...
                n_rel(:, (1:3) + (ii-1)*3) * P_TE((1:3) + (ii-1)*3);
            end
          else
            error('Bad number of n_rel values');
          end

          % Calculate fields in Cartesian coordinates
          E_TE = F_TE * P_TE;
          E_TM = F_TM * P_TM;

          % Apply cartesian to spherical transformation
          % This could also be pre-computed and applied to F
          % but this would be O(Ndipoles) compared to O(2*Nmax),
          % assuming no rotational or mirror symmetry.
          Es_TE = ott.TmatrixDda.apply_cart2sph(E_TE, cart2sph);
          Es_TM = ott.TmatrixDda.apply_cart2sph(E_TM, cart2sph);

          ci = combined_index(n,m);

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

              pq1 = MN(:, modes) \ Es_TE;
              pq2 = MN(:, ~modes) \ Es_TM;

              data(modes,ci) = pq1;
              data(~modes,ci+total_orders) = pq2;

            else

              pq1 = MN\Es_TE;
              pq2 = MN\Es_TM;

              data(:,ci) = pq1;
              data(:,ci+total_orders) = pq2;

            end

          elseif z_rotation > 1
            % Discrete rotational symmetry

            % Calculate which modes preseve symmetry, m = +/- ip
            [alln, allm] = ott.utils.combined_index((1:total_orders).');
            axial_modes = mod(allm - m, z_rotation) == 0;
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

              pq1 = MN(:, modes_evn) \ Es_TE;
              pq2 = MN(:, modes_odd) \ Es_TM;

              data(modes_evn,ci) = pq1;
              data(modes_odd,ci+total_orders) = pq2;

            else

              pq1 = MN(:, modes) \ Es_TE;
              pq2 = MN(:, modes) \ Es_TM;

              data(modes,ci) = pq1;
              data(modes,ci+total_orders) = pq2;

            end

          elseif z_rotation == 0
            % Infinite rotational symmetry
            error('Not yet implemented');
          else
            error('Invalid z_rotation value');
          end

        end		% for n = 1:m
      end		% for m = mrange

    end

    function calc_farfield(tmatrix)
      % TODO: Far-field implementation of T-matrix calculation
      % Warning: Function usage/definition may change without notice
      error('Not yet implemented');
    end
  end
end
