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

      % Calculate voxel locations
      voxels = shape.voxels(spacing);

      % Calculate the T-matrix using DDA
      tmatrix = TmatrixDda(voxels, varargin{:}, ...
        'spacing', spacing);
    end
    
    function k_particle = filter_k_particle(k_particle, mask)
      
      if isscalar(k_particle)
        % Nothing to do
      elseif numel(k_particle) == 3
        % Nothing to do
      elseif numel(k_particle) == numel(mask)
        k_particle(mask) = [];
      elseif numel(k_particle) == 3*numel(mask)
        k_particle(:, mask) = [];
      elseif numel(k_particle) == 9*numel(mask)
        mask = repelem(mask, [3, 3]);
        k_particle(mask) = [];
      else
        error('Unknown filter option');
      end
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
      
      % Check we can allocate sufficient memory
      uV = memory;
      if uV.MaxPossibleArrayBytes < (size(xyz, 2)*3)^2*8
        error('OTT:TmatrixDda:too_many_dipoles', ...
          ['May have too many voxels for calculation, ', ...
          'consider reducing particle size or voxel spacing']);
      end

      % Store inputs k_medium and k_particle
      [k_medium, k_particle] = tmatrix.parser_wavenumber(pa, 2*pi);
      
      % Filter xyz and k_particle for symmetries
      if pa.Results.z_mirror_symmetry
        k_particle = tmatrix.filter_k_particle(k_particle, xyz(3, :) < 0);
        xyz(:, xyz(3, :) < 0) = [];
      end
      if pa.Results.z_rotational_symmetry == 0
        error('z_rotational_symmetry == 0 not yet supported');
      elseif pa.Results.z_rotational_symmetry == 1
        % Nothing to do
      elseif pa.Results.z_rotational_symmetry == 2
        k_particle = tmatrix.filter_k_particle(k_particle, xyz(2, :) < 0);
        xyz(:, xyz(2, :) < 0) = [];
      elseif pa.Results.z_rotational_symmetry == 4
        k_particle = tmatrix.filter_k_particle(...
            k_particle, xyz(2, :) < 0 | xyz(1, :) < 0);
        xyz(:, xyz(2, :) < 0 | xyz(1, :) < 0) = [];
      else
        error('Only z_rotational_symmetry == 2|4 supported for now');
      end

      rtp = ott.utils.xyz2rtp(xyz.');

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
        alpha = pa.Results.polarizability;
      else
        if isempty(pa.Results.spacing)
          error('Spacing is needed for polarizability calculation');
        end
        
        % Get spacing in units of medium wavelength
        spacing = pa.Results.spacing;
        spacing = spacing / wavelength_medium;
        
        switch pa.Results.polarizability
          case 'LDR'
            alpha = ott.utils.polarizability_LDR(...
                spacing, n_relative);
          case 'FCD'
            alpha = ott.utils.polarizability_FCD(...
                spacing, n_relative);
          case 'CM'
            alpha = ott.utils.polarizability_CM(...
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

    function M = cart2sph_mat(theta, phi)

      M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
           cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
           -sin(phi) cos(phi) 0];
         
    end

    function M = sph2cart_mat(theta, phi)

      M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
           cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
           -sin(phi) cos(phi) 0]';
    end
    
    function F = nearfield_matrix_combined(F_total, ...
          m, z_mirror, z_rotation)
      % Combine layers of nearfield_matrix
      %
      % Warning: Function usage/definition may change without notice
      
      F = F_total(:, :, 1);
      
      for ii = 2:z_rotation
        F = F_total(:, :, ii) .* exp(1i*m*2*pi*(ii-1)/z_rotation);
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
            M = eye(3);
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

          % reformat into one column
          M1nm = ott.utils.col3to1(M1);
          N1nm = ott.utils.col3to1(N1);

          MN(:,ci) = -M1nm;
          MN(:,ci+total_orders) = -N1nm;
        end
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
      nrot = 0;

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

      A = ott.utils.interaction_A(k, rtp, alpha);
      
      assert(all(isfinite(A(:))), 'singular polarizability not yet supported');

      F_total = ott.TmatrixDda.nearfield_matrix(...
          theta,phi,r_near, k, xyz, rtp, z_mirror, z_rotation);

      MN = ott.TmatrixDda.calculate_nearfield_modes(...
        k*r_near, theta, phi, Nmax);

      % Allocate memory for T-matrix
      data = zeros(2*total_orders, 2*total_orders);

			for m = mrange
        
        progress_callback(struct('m', m, 'mrange', mrange));
        
        F = ott.TmatrixDda.nearfield_matrix_combined(F_total, ...
            m, z_mirror, z_rotation);
        
        %% TODO: rotsym_interaction_A
        
				for n = max(1, abs(m)):Nmax % step through n incident
          
          [Ei_TE, Ei_TM] = ott.utils.E_inc_vswf(n,-m,rtp,k);

					if iterative
						[P_TE, ~] = gmres(A,Ei_TE,1,1e-6,30);
						[P_TM, ~] = gmres(A,Ei_TM,1,1e-6,30);
					else
						P_TE = A\Ei_TE;
						P_TM = A\Ei_TM;
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
          
          %% TODO: Quad to full

					E_TE = ott.utils.col1to3(F*P_TE);
					E_TM = ott.utils.col1to3(F*P_TM);

					Es_TE2 = zeros(size(E_TE));
					Es_TM2 = zeros(size(E_TM));

					% convert to spherical coordinate vectors
					for j = 1:length(E_TE)
						M = [sin(theta(j))*cos(phi(j)) sin(theta(j))*sin(phi(j)) cos(theta(j));
							cos(theta(j))*cos(phi(j)) cos(theta(j))*sin(phi(j)) -sin(theta(j));
							-sin(phi(j)) cos(phi(j)) 0];

						E = M*E_TE(j,:).';
						Es_TE2(j,:) = E.';

						E = M*E_TM(j,:).';
						Es_TM2(j,:) = E.';
          end
          
          %% TODO: Quad to full again

					% format them to 1 column
					Es_TE = ott.utils.col3to1(Es_TE2);
					Es_TM = ott.utils.col3to1(Es_TM2);

					ci = combined_index(n,m);
          
          % TODO: Check the discrete rotational symettry stuff (nrot)

					if nrot < 2
						pq1 = MN\Es_TE;
						pq2 = MN\Es_TM;
						if ott.utils.isodd(m)
							% shouldn't have to do this but this is to correct a sign error
							data(:,ci) = -pq1;
							data(:,ci+total_orders) = -pq2;
						else
							data(:,ci) = pq1;
							data(:,ci+total_orders) = pq2;
						end
					else
						[~, ~, ci1] = nm_rot(n,m,Nmax,nrot);
						[~, ~, ci2] = nm_rot(n+1,m,Nmax,nrot);
						pq1 = [MN(:,ci1) MN(:,ci2+total_orders)]\Es_TE;
						pq2 = [MN(:,ci2) MN(:,ci1+total_orders)]\Es_TM;
						if ott.utils.isodd(m)
							% shouldn't have to do this but this is to correct a sign error
							data(:,ci) = -pq_col(pq1,ci1,ci2,total_orders);
							data(:,ci+total_orders) = -pq_col(pq2,ci2,ci1,total_orders);
						else
							data(:,ci) = pq_col(pq1,ci1,ci2,total_orders);
							data(:,ci+total_orders) = pq_col(pq2,ci2,ci1,total_orders);
						end
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
