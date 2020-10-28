classdef Dda < ott.tmatrix.Tmatrix
% Construct a T-matrix using the discrete dipole approximation.
% Inherits from :class:`ott.tmatrix.Tmatrix`.
%
% This method calculates how each VSWF mode is scattered and then uses
% point-matching to calculate each column of the T-matrix.
% The field calculations are done using the discrete dipole approximation,
% but the same approach could be used for T-matrix calculation via any
% other field calculation method.
%
% The current DDA implementation requires a lot of memory.  Most small
% desktop computers will be unable to calculate T-matrices for large
% particles (i.e., particles larger than a couple of wavelengths in
% diameter using 20 dipoles per wavelength).  The aim of the next OTT
% release is to include geometric optics and finite difference time domain
% as alternative methods for force calculation with these larger particles.
%
% Properties
%   - dda      -- DDA instance used for field calculation
%   - pmrtp    -- (3xN numeric) Locations for point matching
%
% Static methods
%   - FromShape    -- Construct from a geometric shape
%   - DefaultPmrtp -- Build default ``pmrtp`` locatiosn for point matching
%   - DefaultProgressCallback -- Default progress callback for method
%
% See also :meth:`Dda`, :pkg:`ott.tmatrix.dda`.

  properties (SetAccess=protected)
    dda       % DDA instance used for field calculation
    pmrtp     % (3xN numeric) Locations for point matching
  end

  methods (Static)
    function DefaultProgressCallback(data)
      % Default progress callback for Dda
      %
      % Prints the progress to the terminal.
      %
      % Usage
      %   DefaultProgressCallback(data)
      %
      % Parameters
      %   - data (struct) -- Structure with two fields: index and total.

      disp(['Iteration ' num2str(data.index) ' / ' num2str(max(data.total))]);
    end

    function varargout = FromShape(shape, varargin)
      % Construct a T-matrix from a geometric shape
      %
      % Usage
      %   tmatrix = ott.tmatrix.Dda.FromShape(shape, ...)
      %
      %   [tmatrix, incData, pmData] = ott.tmatrix.Dda.FromShape(...)
      %   Returns the field calculation data for repeated calculations.
      %
      % Optional named arguments
      %   - spacing -- (numeric) -- Dipole spacing in wavelength units.
      %     Default: ``1/20``.
      %
      %   - polarizability -- (function_handle | 3x3 numeric) Method to
      %     calculate polarizability or 3x3 tensor for homogeneous material.
      %     Default: ``@(xyz, spacing, ri) polarizability.LDR(spacing, ri)``
      %
      %   - relative_index -- (function_handle | numeric) Method to calculate
      %     relative refractive index or homogeneous value.  Ignored if
      %     polarizability is a 3x3 tensor.
      %
      %   - low_memory -- (logical) If we should use the low memory DDA
      %     implementation.  Default: ``false``.
      %
      % Additional parameters passed to class constructor.

      p = inputParser;
      p.addParameter('spacing', 1/20, @isnumeric);
      p.addParameter('polarizability', ...
          @(~, s, r) ott.tmatrix.dda.polarizability.LDR(s, r));
      p.addParameter('relative_index', 1.0);
      p.addParameter('low_memory', false);
      p.addParameter('xySymmetry', shape.xySymmetry);
      p.addParameter('zRotSymmetry', shape.zRotSymmetry);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct a DDA representation
      dda = ott.tmatrix.dda.Dda.FromShape(shape, ...
          'spacing', p.Results.spacing, ...
          'polarizability', p.Results.polarizability, ...
          'relative_index', p.Results.relative_index, ...
          'xySymmetry', p.Results.xySymmetry, ...
          'zRotSymmetry', p.Results.zRotSymmetry);

      % Select the high memory implementation
      if ~p.Results.low_memory
        dda = ott.tmatrix.dda.DdaHighMem(dda);
      end

      % Calculate the T-matrix
      [varargout{1:nargout}] = ott.tmatrix.Dda(dda, unmatched{:});
    end

    function pmrtp = DefaultPmrtp(Nmax, varargin)
      % Build default grid of points for point matching
      %
      % Usage
      %   pmrtp = DefaultPmrtp(Nmax, ...)
      %
      % Parameters
      %   - Nmax -- (numeric) Row Nmax for generated T-matrix.
      %
      % Optional named parameters
      %   - radius -- (numeric) Radius for point matching locations.
      %     Must be positive scalar or Inf for far-field point matching.
      %     Default: ``Inf``.
      %
      %   - angulargrid -- ({theta, phi}) Angular grid of points for
      %     calculation of radii.  Default is equally spaced angles with
      %     the number of points determined by Nmax.
      %
      %   - xySymmetry -- (logical) If the generated grid should be for
      %     mirror symmetry DDA.  Default: ``false``.
      %
      %   - zRotSymmetry -- (numeric) If the generated grid should be for
      %     rotationally symmetric DDA.  Default: ``1``.

      p = inputParser;
      p.addParameter('radius', Inf, @isnumeric);
      p.addParameter('angulargrid', []);
      p.addParameter('xySymmetry', false);
      p.addParameter('zRotSymmetry', 1);
      p.parse(varargin{:});

      % Get or calculate angular grid
      angulargrid = p.Results.angulargrid;
      if isempty(angulargrid)

        % We need at least one point for every beam shape coefficient
        % With no mirror or rotational symmetry, use the same grid as PmGauss
        ntheta = (Nmax + 1);
        nphi = 2*(Nmax + 1);

        % We can reduce the number of point around the z-axis when we
        % have z rotational symmetry since modes will only scatter to
        % other modes with a multiple of the rotational symmetry factor.
        if p.Results.zRotSymmetry > 1
          nphi = ceil(nphi ./ p.Results.zRotSymmetry);
        elseif p.Results.zRotSymmetry == 0

          % TODO: When we work out how to infinite rotational DDA
          %   we should change this value to 1, for now its 3.
          nphi = 3;
        end

        % Similarly, we can reduce the number of point in theta
        % when we have z-mirror symmetry since we match twice
        if p.Results.xySymmetry
          ntheta = ceil(ntheta ./ 2);
        end

        theta = ((1:ntheta)-0.5)/ntheta * pi;
        phi = ((1:nphi)-1)/nphi * 2*pi;

        [theta, phi] = meshgrid(theta, phi);

        % We also need to rescale our points by a similar amount
        % This avoids making the problem rank deficient
        if p.Results.zRotSymmetry > 1
          phi = phi ./ p.Results.zRotSymmetry;
        end
        if p.Results.xySymmetry
          theta = theta ./ 2;
        end

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

      % Create array of angular grid points
      pmrtp = [p.Results.radius*ones(size(theta(:))), theta(:), phi(:)].';
    end
  end

  methods
    function [tmatrix, incData, pmData] = Dda(dda, varargin)
      % Construct a T-matrix using a DDA simulation for field calculation
      %
      % Usage
      %   tmatrix = Dda(dda, ...)
      %
      %   [tmatrix, incData, pmData] = Dda(dda, ...)
      %   Also returns the VSWF data structures used for calculating the
      %   incident field and the point matching field.
      %
      % Parameters
      %   - dda -- (ott.tmatrix.dda.Dda instance) The DDA instance
      %     used for field calculations.
      %
      % Optional named arguments
      %   - Nmax -- (numeric) Size of the VSWF expansion used for the
      %     T-matrix point matching (determines T-matrix rows).
      %     Default: ``ott.utils.ka2nmax(2*pi*shape.maxRadius)`` (may
      %     need different values to give convergence for some shapes).
      %
      %   - ci -- (N numeric) Number of modes to calculate scattering for.
      %     This determines number of T-matrix columns.
      %     Default: ``1:ott.utils.combined_index(Nmax, Nmax)`` (all modes).
      %
      %   - pmrtp -- (3xN numeric) Coordinatse for point matching.
      %     Radial coordinate must all be finite or all be Inf.
      %     Default uses ``ott.tmatrix.Dda.DefaultPmrtp``.
      %
      %   - incData (ott.utils.VswfData) -- Data structure for repeated
      %     incident field calculations.
      %     Default: ``ott.utils.VswfData()``.
      %
      %   - pmData (ott.utils.VswfData) -- Data structure for repeated
      %     point matching field calculations.
      %     Default: ``ott.utils.VswfData()``.
      %
      % Unmatched parameters passed to :meth:`calculate_columns`.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('Nmax', []);
      p.addParameter('ci', []);
      p.addParameter('pmrtp', []);
      p.addParameter('incData', ott.utils.VswfData());
      p.addParameter('pmData', ott.utils.VswfData());
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Store dda instance
      tmatrix.dda = dda;

      % Get or calculate Nmax
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        rmax = max(vecnorm(dda.locations));
        Nmax = ott.utils.ka2nmax(2*pi*rmax);
      else
        assert(isnumeric(Nmax) && isscalar(Nmax) && Nmax > 0, ...
            'Nmax must be positive numeric scalar');
      end

      % Get or calculate ci for calculation
      ci = p.Results.ci;
      if isempty(ci)
        ci = 1:ott.utils.combined_index(Nmax, Nmax);
      else
        assert(isnumeric(ci) && isvector(ci) && all(ci > 0), ...
            'ci must be positive numeric vector');
      end

      % Get or calculate pmrtp
      pmrtp = p.Results.pmrtp;
      if isempty(pmrtp)
        pmrtp = ott.tmatrix.Dda.DefaultPmrtp(Nmax, ...
            'xySymmetry', dda.xySymmetry, ...
            'zRotSymmetry', dda.zRotSymmetry);
      else
        assert(ismatrix(pmrtp) && isnumeric(pmrtp) ...
            && size(pmrtp, 1) == 3, ...
            'pmrtp must be 3xN numeric matrix');
        assert(all(pmrtp(1, :) == pmrtp(1, 1)) && pmrtp(1, 1) > 0, ...
            'prtp(1, :) must be constant positive value');
      end
      tmatrix.pmrtp = pmrtp;

      % Calculate column Nmax
      [nmodes, ~] = ott.utils.combined_index(ci);
      cNmax = max(nmodes);

      % Calculate total number of orders
      total_orders_cols = ott.utils.combined_index(cNmax, cNmax);
      total_orders_rows = ott.utils.combined_index(Nmax, Nmax);

      % Allocate memory for T-matrix data
      tmatrix.data = sparse(2*total_orders_rows, 2*total_orders_cols);
      tmatrix = tmatrix.setType('scattered');

      % Calculate columns
      [tmatrix, incData, pmData] = tmatrix.calculate_columns(ci, ...
          p.Results.incData, p.Results.pmData, unmatched{:});
    end

    function [tmatrix, incData, pmData] = calculate_columns(tmatrix, ...
        ci, incData, pmData, varargin)
      % Calculate T-matrix columns
      %
      % This function is called by the class constructor but could also
      % be called separately for parallel T-matrix calculation.
      %
      % Usage
      %   [tmatrix, incData, pmData] = tmatrix.calculate_columns(...
      %       ci, incData, pmData, ...)
      %
      % Parameters
      %   - ci -- (N numeric) Combined indices for columns.
      %
      %   - incData -- (ott.utils.VswfData) Data structure for VSWF
      %     incident field calculations.
      %
      %   - pmData -- (ott.utils.VswfData) Data structure for point
      %     matching VSWF field calculations.
      %
      % Optional named parameters
      %   - solver -- (function_handle) Solver to use.  Good solvers to try
      %     include ``gmres``, ``bicgstab`` and ``\``.
      %     Default: ``@(A, E) A \ E``.
      %
      %   - low_memory -- (logical) If the method should use the reduced
      %     memory dipole field calculation.
      %     Default: ``false`` if ``dda`` is a :class:`DdaHighMem` instance
      %     and ``true`` otherwise.
      %
      %   - progress (function_handle) -- Function to call for progress
      %     updates during method evaluation.  Takes one argument, see
      %     :meth:`DefaultProgressCallback` for more information.
      %     Default: ``@DefaultProgressCallback``.

      p = inputParser;
      p.addParameter('low_memory', ...
          ~isa(tmatrix.dda, 'ott.tmatrix.dda.DdaHighMem'));
      p.addParameter('progress', @ott.tmatrix.Dda.DefaultProgressCallback);
      p.addParameter('solver', @(A, E) A \ E, @(x) isa(x, 'function_handle'));
      p.parse(varargin{:});

      ott.utils.nargoutCheck(tmatrix, nargout);

      % Pre-calculate incident fields
      [Ei, incData] = tmatrix.calculate_incident_field(ci, incData);

      % Pre-calculate outgoing fields for PM
      [MN, pmData] = tmatrix.calculate_outgoing_field(pmData);

      % Pre-calculate transfer matrix for faster calculations
      F = tmatrix.calculate_scattered_matrix(p.Results.low_memory);

      [nmodes, mmodes] = ott.utils.combined_index(ci(:).');

      for m = unique(mmodes)

        % Report progress
        p.Results.progress(struct(...
          'index', sum(mmodes < m), 'total', numel(mmodes)));

        % Calculate which nmodes/Ei we need for calculation
        ournmodes = nmodes(mmodes == m);
        ourEi = Ei(:, [mmodes == m, mmodes == m]);

        % Calculate scattered fields
        En = tmatrix.calculate_scattered_fields(m, ournmodes, ...
          ourEi, F, p.Results.solver);

        % Do point matching for each column
        tmatrix = tmatrix.calculate_pm_columns(m, ournmodes, MN, En);
      end

      % Final Report progress
      p.Results.progress(struct(...
        'index', numel(mmodes), 'total', numel(mmodes)));
    end
  end

  methods (Hidden)
    function [Ei, data] = calculate_incident_field(tmatrix, ci, data)
      % Pre-calculate the incident field data
      %
      % Usage
      %   [Ei, data] = tmatrix.calculate_incident_field(ci, data)

      % Calculate VSWF basis set
      beam = ott.bsc.Bsc.BasisSet(ci);

      % Calculate FieldVector data
      [fE, data] = beam.efield(tmatrix.dda.locations.*1.2, ...
          'data', data, 'basis', 'regular');

      % Package for DDA
      Ei = reshape(fE.vxyz, 3*size(tmatrix.dda.locations, 2), 2*numel(ci));
    end

    function [MN, data] = calculate_outgoing_field(tmatrix, data)
      % Pre-calculate the outgoing fields and transfer matrix for PM
      %
      % Usage
      %   [MN, data] = tmatrix.calculate_outgoing_field(data)

      % Get modes to match
      ci = 1:size(tmatrix.data, 1)/2;

      % Calculate VSWF basis set
      beam = ott.bsc.Bsc.BasisSet(ci);

      % Calculate Field vectors
      if isfinite(tmatrix.pmrtp(1, 1))
        [fMN, data] = beam.efieldRtp(tmatrix.pmrtp, 'data', data, ...
            'basis', 'outgoing');
        fMN = fMN.vxyz.*2;
      else
        [fMN, data] = beam.efarfield(tmatrix.pmrtp, 'data', data, ...
            'basis', 'outgoing');
        fMN = fMN.vxyz./(2*pi);
      end

      % Package
      MN = reshape(fMN, 3*size(tmatrix.pmrtp, 2), 2*numel(ci));
    end

    function F = calculate_scattered_matrix(tmatrix, low_memory)
      % Pre-calculate dipole scattered field matrix if not low_memory

      if low_memory
        F = [];  % Nothing to do
      else
        dipoles = ott.tmatrix.dda.Dipole(tmatrix.dda.locations, ...
          zeros(3*size(tmatrix.dda.locations, 2), 0), ...
          'xySymmetry', tmatrix.dda.xySymmetry, ...
          'zRotSymmetry', tmatrix.dda.zRotSymmetry);

        if isfinite(tmatrix.pmrtp(1, 1))
          xyz = ott.utils.rtp2xyz(tmatrix.pmrtp);
          F = dipoles.enearfield_matrix(xyz, 'low_memory', false);
        else
          F = dipoles.efarfield_matrix(tmatrix.pmrtp, 'low_memory', false);
        end
      end
    end

    function E = calculate_scattered_fields_internal(tmatrix, dipoles, F)
      % Calculate scattered fields at specific locations

      if isempty(F)
        if isfinite(tmatrix.pmrtp(1, 1))
          xyz = ott.utils.rtp2xyz(tmatrix.pmrtp);
          E = dipoles.efield(xyz, 'low_memory', true);
        else
          E = dipoles.efarfield(tmatrix.pmrtp, 'low_memory', true);
        end
        
        % Reshape to match F*dipoles output
        E = reshape(E.vxyz, 3*size(tmatrix.pmrtp, 2), []);
      else
        E = F * dipoles;
      end
    end

    function En = calculate_scattered_fields(tmatrix, ...
        m, nmodes, Ei, F, solver)
      % Solve DDA problem and calculate scattered fields

      if tmatrix.dda.xySymmetry
        evn = mod(m + nmodes, 2) == 0;

        dipolesOdd = tmatrix.dda.solve(Ei(:, [evn, ~evn]), ...
            'parity', 'odd', 'rorder', m, 'solver', solver);
        dipolesEvn = tmatrix.dda.solve(Ei(:, [~evn, evn]), ...
            'parity', 'even', 'rorder', m, 'solver', solver);

        Eodd = tmatrix.calculate_scattered_fields_internal(dipolesOdd, F);
        Eevn = tmatrix.calculate_scattered_fields_internal(dipolesEvn, F);

        En = zeros(size(Eodd, 1), 2*numel(nmodes));
        En(:, [evn, evn]) = [Eodd(:, 1:sum(evn)), Eevn(:, sum(~evn)+1:end)];
        En(:, ~[evn, evn]) = [Eevn(:, 1:sum(~evn)), Eodd(:, sum(evn)+1:end)];

      else
        dipoles = tmatrix.dda.solve(Ei, 'rorder', m, 'solver', solver);
        En = tmatrix.calculate_scattered_fields_internal(dipoles, F);
      end
    end

    function tmatrix = calculate_pm_columns(tmatrix, m, nmodes, MN, En)
      % Calculate individual column pair of T-matrix

      ott.utils.nargoutCheck(tmatrix, nargout);

      total_orders_cols = size(tmatrix.data, 2)/2;
      total_orders_rows = size(tmatrix.data, 1)/2;

      [alln, allm] = ott.utils.combined_index((1:total_orders_rows).');

      % Find even modes
      even_modes = logical(mod(alln + allm, 2));

      z_rotation = tmatrix.dda.zRotSymmetry;
      z_mirror = tmatrix.dda.xySymmetry;

      for n = nmodes

        % Calculate column index
        ci = ott.utils.combined_index(n, m);

        % Get mode fields
        E_TE = En(:, [nmodes == n, false(1, numel(nmodes))]);
        E_TM = En(:, [false(1, numel(nmodes)), nmodes == n]);

        if z_rotation == 1
          % No rotational symmetry

          if z_mirror
            % Only mirror symmetry
            % Parity is conserved, even modes go to even modes, etc.
            % Reference: https://arxiv.org/pdf/physics/0702045.pdf

            modes = [ even_modes; ~even_modes ];

            if ~logical(mod(n + m, 2))
              modes = ~modes;
            end

            pq1 = MN(:, modes) \ E_TE;
            pq2 = MN(:, ~modes) \ E_TM;

            tmatrix.data(modes,ci) = pq1;
            tmatrix.data(~modes,ci+total_orders_cols) = pq2;

          else

            pq1 = MN\E_TE;
            pq2 = MN\E_TM;

            tmatrix.data(:,ci) = pq1;
            tmatrix.data(:,ci+total_orders_cols) = pq2;

          end

        else

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
            if logical(mod(n + m, 2))
              modes_evn = modes & [ even_modes; ~even_modes ];
              modes_odd = modes & [ ~even_modes; even_modes ];
            else
              modes_odd = modes & [ even_modes; ~even_modes ];
              modes_evn = modes & [ ~even_modes; even_modes ];
            end

            pq1 = MN(:, modes_evn) \ E_TE;
            pq2 = MN(:, modes_odd) \ E_TM;

            tmatrix.data(modes_evn,ci) = pq1;
            tmatrix.data(modes_odd,ci+total_orders_cols) = pq2;

          else

            pq1 = MN(:, modes) \ E_TE;
            pq2 = MN(:, modes) \ E_TM;

            tmatrix.data(modes,ci) = pq1;
            tmatrix.data(modes,ci+total_orders_cols) = pq2;

          end
        end

      end
    end
  end
end
