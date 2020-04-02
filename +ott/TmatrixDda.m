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
        'z_mirror_symmetry', z_mirror_symmetry & false);
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
      %
      %   - verbose (logical) -- Display additional information.
      %     Doesn't affect the display of the progress callback.
      %     Default: false.

      import ott.TmatrixDda;
      import ott.Tmatrix;

      tmatrix = tmatrix@ott.Tmatrix();

      % Parse inputs
      pa = TmatrixDda.input_parser(varargin{:});
      
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
        error('Not yet implemented');
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
      if uV.MaxPossibleArrayBytes < (size(xyz, 2)*3)^2*8
        error('OTT:TmatrixDda:too_many_dipoles', ...
          ['May have too many voxels for calculation, ', ...
          'consider reducing particle size or voxel spacing']);
      end

      rtp = ott.utils.xyz2rtp(xyz.');
      
      % Tell the user some things
      if pa.Results.verbose && pa.Results.z_rotational_symmetry ~= 1
        disp(['Voxels with rotationa symmetry: ' num2str(size(xyz, 2))]);
      end

      % Get or estimate Nmax from the inputs
      Nmax = pa.Results.Nmax;
      if isempty(Nmax)
        
        % Get or estimate spacing
        spacing = pa.Results.spacing;
        if isempty(spacing)
          warning('Estimating spacing for Nmax calculation from k_particle');
          spacing = 2*pi./max(abs(k_particle))./20;
        end
        
        maxRadius = max(rtp(:, 1)) + spacing/2;
        Nmax = ott.utils.ka2nmax(maxRadius * abs(k_medium));
      end
      assert(isscalar(Nmax), 'Only scalar Nmax supported for now');

      % Calculate range of m-orders
      mrange = -Nmax:Nmax;

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

      % Calculate the T-matrix
      tmatrix.data = tmatrix.calc_nearfield(Nmax, xyz.', rtp, ...
          n_relative, mrange, alpha, pa.Results.progress_callback, ...
          pa.Results.z_mirror_symmetry, pa.Results.z_rotational_symmetry);

      % Store the type of T-matrix
      tmatrix.type = 'scattered';

    end
  end

  methods (Hidden, Static)
    
    function k_particle = remove_by_mask(k_particle, mask)
      
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

      M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
           cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
           -sin(phi) cos(phi) 0];
         
    end

    function M = sph2cart_mat(theta, phi)

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
    
    function F = nearfield_matrix_combined(F_total, ...
          m, z_mirror, z_rotation)
      % Combine layers of nearfield_matrix
      %
      % Warning: Function usage/definition may change without notice
      
      F = F_total(:, :, 1);
      
      for ii = 2:z_rotation
        F = F + F_total(:, :, ii) .* exp(1i*m*2*pi*(ii-1)/z_rotation);
      end
    end

    function F = nearfield_matrix(theta, phi, r_near, k, xyz, rtp, ...
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
          
          % Calculate coordinate rotation
          if m == 1
            M = 1;
          else
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
    
    function A = interaction_A(k,r,varargin)
      % Calculate the interaction matrix
      %
      % TODO: Should this be part of TmatrixDda?
      %
      % A = interaction_A(k_medium, voxels, alpha) calculate the off-diagonal
      % terms and use the polarisability, alpha, to form the diagonal.
      % voxels must be a Nx3 array of Cartesian voxel coordinates.
      %
      % A = interaction_A(k_medium, voxels, 'inv_alpha', inv_alpha) as above,
      % but uses the inverse polarisability.
      %
      % alpha and inv_alpha can be scalar and 3x1 vectors or 3x3 matrix for
      % homogeneous isotropic and birefringent material.  For inhomogeneous
      % materials, use N, 3xN (or 3N vector) or 3x3N matrices respectively.
      % If alpha or inv_alpha are not supplied, the diagonal is left empty.
      %
      % Based on OMG/Code/DDA/interaction_A

      % Copyright 2018 Isaac Lenton

      p = inputParser;
      p.addOptional('alpha', []);
      p.addParameter('inv_alpha', []);
      p.addParameter('z_rotational_symmetry', 1);
      p.parse(varargin{:});

      if ~isempty(p.Results.alpha) && ~isempty(p.Results.inv_alpha)
        error('Either alpha or inv_alpha must be supplied or none, not both');
      end

      assert(p.Results.z_rotational_symmetry >= 1, ...
        'Only discrete rotational symmetry supported for now');

      assert(ismatrix(r) && size(r, 2) == 3, 'voxels must be a Nx3 matrix');

      N = size(r, 1);

      % Pre-compute cart2sph transformations for each dipole
      cart2sph = [];
      sph2cart = [];
      if p.Results.z_rotational_symmetry > 1

        % Pre-allocate memory
        cart2sph = zeros(3, 3*N);
        sph2cart = cart2sph;

        % Compute
        for ii = 1:N
          [~, theta, phi] = ott.utils.xyz2rtp(r(ii, :));

          cart2sph(:, (1:3) + (ii-1)*3) = ott.TmatrixDda.cart2sph_mat(...
              theta, phi);
        end
      end

      % Generate the matrix of off-diagonal elements
      A = zeros(3*N,3*N, p.Results.z_rotational_symmetry);
      for m = 1:p.Results.z_rotational_symmetry

        % Pre-compute sph2cart transformations for each dipole
        if m > 1
          for ii = 1:N
            [~, theta, phi] = ott.utils.xyz2rtp(r(ii, :));

            sph2cart(:, (1:3) + (ii-1)*3) = ott.TmatrixDda.sph2cart_mat(...
                theta, phi + 2*pi*(m-1)/p.Results.z_rotational_symmetry);
          end
        end

        for j=1:N
          A((1:3) + 3*(j-1),:, m) = ott.TmatrixDda.calc_Aj(k,r,j, p.Results.z_rotational_symmetry, m, ...
            cart2sph, sph2cart);
        end
      end

      % Put the inverse polarisability elements into usable form
      inv_alpha = [];
      if ~isempty(p.Results.alpha)
        val = p.Results.alpha;
        sz = size(val);
        if all(sz == [1, 1])
          inv_alpha = repmat([1./val, 0, 0; 0, 1./val, 0; 0, 0, 1./val], 1, N);
        elseif all(sz == [3, 3])
          inv_alpha = repmat(inv(val), 1, N);
        elseif all(sz == [3, 1])
          inv_alpha = repmat(diag(1.0./val), 1, N);
        elseif numel(val) == N
          inv_alpha = zeros(3, 3*N);
          inv_alpha(1, 1:3:end) = 1.0./val;
          inv_alpha(2, 2:3:end) = 1.0./val;
          inv_alpha(3, 3:3:end) = 1.0./val;
        elseif all(sz == [3, N])
          inv_alpha = zeros(3, 3*N);
          inv_alpha(1, 1:3:end) = 1.0./val(1:3:end);
          inv_alpha(2, 2:3:end) = 1.0./val(2:3:end);
          inv_alpha(3, 3:3:end) = 1.0./val(3:3:end);
        elseif all(sz == [3, 3*N])
          inv_alpha = zeros(3, 3*N);
          for ii = 1:N
            inv_alpha(:, (1:3) + 3*(ii-1)) = inv(val(:, (1:3) + 3*(ii-1)));
          end
        else
          error('Unsupported size of alpha');
        end
      elseif ~isempty(p.Results.inv_alpha)
        val = p.Results.inv_alpha;
        sz = size(val);
        if all(sz == [1, 1])
          inv_alpha = repmat([val, 0, 0; 0, val, 0; 0, 0, val], 1, N);
        elseif all(sz == [3, 3])
          inv_alpha = repmat(val, 1, N);
        elseif all(sz == [3, 1])
          inv_alpha = repmat(diag(val), 1, N);
        elseif numel(val) == N
          inv_alpha = zeros(3, 3*N);
          inv_alpha(1, 1:3:end) = val;
          inv_alpha(2, 2:3:end) = val;
          inv_alpha(3, 3:3:end) = val;
        elseif all(sz == [3, N])
          inv_alpha = zeros(3, 3*N);
          inv_alpha(1, 1:3:end) = val(1:3:end);
          inv_alpha(2, 2:3:end) = val(2:3:end);
          inv_alpha(3, 3:3:end) = val(3:3:end);
        elseif all(sz == [3, 3*N])
          inv_alpha = val;
        else
          error('Unsupported size of inv_alpha');
        end
      end

      % Put the inverse polarisability on the diagonal
      if ~isempty(inv_alpha)
        Ac = mat2cell(inv_alpha, 3, repmat(3,1,N));
        A(:, :, 1) = A(:, :, 1) + blkdiag(Ac{:});
      end

    end

    function Aj = calc_Aj(k_medium, r, j, z_rotation, m, cart2sph, sph2cart)
      % calculates a 3 X 3N block comprising N number of 3 X 3 Green's tensors

      % The following method is about twice as slow as the previous
      % implementation (for no rot. sym.), but it has a lot of similarity
      % with the F matrix and perhaps we can write a nice & optimal version
      % later!!!

      % Get voxel location for this quadrant
      if m == 1
        our_xyz = r;
      else
        rotphi = (m-1) * 2*pi/z_rotation;
        rotM = [cos(rotphi) -sin(rotphi) 0; sin(rotphi) cos(rotphi) 0; 0 0 1];
        our_xyz = (rotM * (r.')).';
      end

      n_dipoles = size(r, 1);

      r_vec = our_xyz - r(j, :);
      r_jk = vecnorm(r_vec, 2, 2);

      % Remove diagonal term, only applies when not calculating
      % elements with rotational symmetry.  (avoids nan)
      if m == 1
        r_jk(j) = 1;
      end

      r_hat = r_vec./r_jk;

      Aj = zeros(3, 3*n_dipoles);

      for ii = 1:n_dipoles

        % Skip self-interaction term
        if ii == j && m ~= 1
          continue;
        end

        % Calculate coordinate rotation
        if m == 1
          M = 1;
        else
          M = sph2cart(:, (1:3) + (ii-1)*3) * cart2sph(:, (1:3) + (ii-1)*3);
        end

        rr = r_hat(ii, :).'*r_hat(ii, :);

        Aj(:, (1:3) + 3*(ii-1)) = exp(1i*k_medium*r_jk(ii, :))/r_jk(ii, :)*...
          (k_medium^2*(rr - eye(3)) + (1i*k_medium*r_jk(ii, :)-1)/r_jk(ii, :)^2*(3*rr - eye(3))) * M;

      end
    end

    function data = calc_nearfield(Nmax, xyz, rtp, n_rel, mrange, alpha, ...
        progress_callback, z_mirror, z_rotation)
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

      A_total = ott.TmatrixDda.interaction_A(k, rtp, alpha, ...
        'z_rotational_symmetry', z_rotation);
      
      assert(all(isfinite(A_total(:))), 'singular polarizability not yet supported');

      F_total = ott.TmatrixDda.nearfield_matrix(...
          theta,phi,r_near, k, xyz, rtp, z_mirror, z_rotation);

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

			for m = mrange
        
        progress_callback(struct('m', m, 'mrange', mrange));
        
        F = ott.TmatrixDda.nearfield_matrix_combined(F_total, ...
            m, z_mirror, z_rotation);
        
        A = ott.TmatrixDda.nearfield_matrix_combined(A_total, ...
            m, z_mirror, z_rotation);
        
				for n = max(1, abs(m)):Nmax % step through n incident
          
          [Ei_TE, Ei_TM] = ott.utils.vswfcart(n, m, ...
              rtp(:, 1)*k, rtp(:, 2), rtp(:, 3), 'regular');
          Ei_TE = Ei_TE.';
          Ei_TM = Ei_TM.';

          if iterative
            [P_TE, ~] = gmres(A,Ei_TE(:),1,1e-6,30);
            [P_TM, ~] = gmres(A,Ei_TM(:),1,1e-6,30);
          else
            P_TE = A\Ei_TE(:);
            P_TM = A\Ei_TM(:);
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
          E_TE = F * P_TE;
          E_TM = F * P_TM;
          
          % Apply cartesian to spherical transformation
          % This could also be pre-computed and applied to F
          % but this would be O(Ndipoles) compared to O(2*Nmax),
          % assuming no rotational or mirror symmetry.
          Es_TE = ott.TmatrixDda.apply_cart2sph(E_TE, cart2sph);
          Es_TM = ott.TmatrixDda.apply_cart2sph(E_TM, cart2sph);

					ci = combined_index(n,m);

					if z_rotation == 1
            % No rotational symmetry
            
						pq1 = MN\Es_TE;
						pq2 = MN\Es_TM;
            
            data(:,ci) = pq1;
            data(:,ci+total_orders) = pq2;
            
          elseif z_rotation > 1
            % Discrete rotational symmetry
            
            % Calculate which modes preseve symmetry, m = +/- ip
            [~, allm] = ott.utils.combined_index((1:total_orders).');
            axial_modes = mod(allm - m, z_rotation) == 0;
            modes = [axial_modes; axial_modes];
            
            pq1 = MN(:, modes) \ Es_TE;
            pq2 = MN(:, modes) \ Es_TM;
            
            data(modes,ci) = pq1;
            data(modes,ci+total_orders) = pq2;
            
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
