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
      % See :meth:`TmatrixDda` for optional named parameters.

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
      tmatrix = TmatrixDda(voxels, varargin{:});
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
      %
      %   - index_relative (numeric)  -- Relative refractive index.
      %     Default: 1.0
      %   - wavelength0 (numeric)     -- Wavelength in vacuum.
      %     Default: 1.0
      %   - spacing (numeric)         -- Spacing between dipoles.
      %     Default: ``wavelength_particle/10``
      %   - polarizability (enum|numeric) -- Polarizability or method
      %     name to use to calculate from relative refractive index.
      %     Default: 'LDR'.  Supported methods: 'LDR', 'FCD', 'CM'.

      import ott.TmatrixDda;
      import ott.Tmatrix;

      tmatrix = tmatrix@ott.Tmatrix();

      % Parse inputs
      pa = TmatrixDda.input_parser(varargin{:});

      % Store inputs k_medium and k_particle
      [k_medium, k_particle] = tmatrix.parser_wavenumber(pa, 2*pi);

      rtp = ott.utils.xyz2rtp(xyz.');

      % Get spacing
      spacing = pa.Results.spacing;
      if isempty(spacing)
        spacing = 2*pi./max(abs(k_particle))./20;
      end

      import ott.utils.ka2nmax;

      % Get or estimate Nmax from the inputs
      Nmax = pa.Results.Nmax;
      if isempty(Nmax)
        maxRadius = max(rtp(:, 1)) + spacing/2;
        Nmax = ka2nmax(maxRadius * abs(k_medium));
      end

      % Calculate range of m-orders
      mrange = -Nmax:Nmax;

      % Put everything in units of wavelengths
      wavelength_medium = 2*pi./k_medium;
      spacing = spacing / wavelength_medium;
      xyz = xyz ./ wavelength_medium;
      rtp(:, 1) = rtp(:, 1) ./ wavelength_medium;

      n_relative = k_particle./k_medium;

      % Compute or get polarizability from inputs
      if isnumeric(pa.Results.polarizability)
        alpha = pa.Results.polarizability;
      else
        switch pa.Results.polarizability
          case 'LDR'
            alpha = ott.utils.polarizability_LDR(...
                spacing, n_relative, 2*pi);
          case 'FCD'
            alpha = ott.utils.polarizability_FCD(...
                spacing, n_relative, 2*pi);
          case 'CM'
            alpha = ott.utils.polarizability_CM(...
                spacing, n_relative);
          otherwise
            error('Unknown polarizability method name');
        end
      end

      % Calculate the T-matrix
      tmatrix.data = tmatrix.calc_nearfield(Nmax, xyz.', rtp, ...
          n_relative, spacing, mrange, alpha, pa.Results.progress_callback);

      % Store the type of T-matrix
      tmatrix.type = 'scattered';

    end
  end

  methods (Access=protected, Static)

    function F = nearfield_matrix(theta, phi, r_near, k, xyz, rtp)
      % Based on DDA/T-matrix/near_field/nearfield_matrix

      N = size(rtp, 1);

      np = length(theta);
      r_E = zeros(np,3);
      [r_E(:,1), r_E(:,2), r_E(:,3)] = ott.utils.rtp2xyz(r_near,theta,phi);

      F = zeros(3*np,3*N);
      I = eye(3);
      for p = 1:np

        p_i = 3*(p-1)+1;
        for n = 1:N
          r_vec = (r_E(p,:)-xyz(n,:));
          r_jk = norm(r_vec);
          r_hat = r_vec/r_jk;
          rr = r_hat'*r_hat;

          n_i = 3*(n-1)+1;
          F(p_i:p_i+2,n_i:n_i+2) = exp(1i*k*r_jk)/r_jk*...
            (k^2*(rr - I) + (1i*k*r_jk-1)/r_jk^2*(3*rr - I));
        end
      end
    end

    function data = calc_nearfield(Nmax, xyz, rtp, n_rel, d, mrange, alpha, ...
        progress_callback)
      % Near-field implementation of T-matrix calculation
      %
      % This is based on DDA/T-matrix/near_field/DDA_T_NF.m

      iterative = true;
      k = 2*pi;
      nrot = 0;

      import ott.utils.combined_index;

      total_orders = combined_index(Nmax,Nmax);

      % What is this?
      npts = round(sqrt(8*total_orders));
      [theta,phi] = ott.utils.angulargrid(npts,npts);
      npts = length(theta);

      % Hmm, what if we have a larger sphere?
      r_near = 8; % near field radius

      A = ott.utils.interaction_A(k, rtp, alpha);

      F = ott.TmatrixDda.nearfield_matrix(...
          theta,phi,r_near, k, xyz, rtp);

      % precalculate M1 * N1's for all pts and all modes
      MN = zeros(3*npts,2*total_orders);
      for n = 1:Nmax
        for m = -n:n
          ci = combined_index(n,m);
          [M1,N1] = ott.utils.vswf(n,m,k*r_near,theta,phi,1);

          % reformat into one column
          M1nm = ott.utils.col3to1(M1);
          N1nm = ott.utils.col3to1(N1);

          MN(:,ci) = -M1nm ./ n_rel;
          MN(:,ci+total_orders) = -N1nm;
        end
      end

      % Allocate memory for T-matrix
      data = zeros(2*total_orders, 2*total_orders);

			for m = mrange
        
        progress_callback(struct('m', m, 'mrange', mrange));
        
				for n = max(1, abs(m)):Nmax % step through n incident

% 					Nn = 1/sqrt(n*(n+1));

          [Ei_TE, Ei_TM] = ott.utils.E_inc_vswf(n,-m,rtp,k);

					if iterative
						[P_TE, ~] = gmres(A,Ei_TE,1,1e-6,30);
						[P_TM, ~] = gmres(A,Ei_TM,1,1e-6,30);
					else
						P_TE = A\Ei_TE;
						P_TM = A\Ei_TM;
					end

					E_TE = ott.utils.col1to3(F*P_TE);
					E_TM = ott.utils.col1to3(F*P_TM);

					Es_TE2 = zeros(size(E_TE));
					Es_TM2 = zeros(size(E_TM));

					% convert to spherical coordinate vectors
					for j = 1:length(E_TE)
						M = [sin(theta(j))*cos(phi(j)) sin(theta(j))*sin(phi(j)) cos(theta(j));
							cos(theta(j))*cos(phi(j)) cos(theta(j))*sin(phi(j)) -sin(theta(j));
							-sin(phi(j)) cos(phi(j)) 0];

						E = M*E_TE(j,:)';
						Es_TE2(j,:) = E';

						E = M*E_TM(j,:)';
						Es_TM2(j,:) = E';
					end

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
      error('Not yet implemented');
    end
  end
end
