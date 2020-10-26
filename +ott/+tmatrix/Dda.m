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
%   - pmradius -- Radius for point matching (either Inf or a finite value)
%
% Static methods
%   - FromShape   -- Construct from a geometric shape
%
% See also :meth:`Dda`, :pkg:`ott.tmatrix.dda`.

  properties (SetAccess=protected)
    dda       % DDA instance used for field calculation
    pmradius  % Radius for point matching (either Inf or a finite value)
  end

  methods (Static)
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
      %   - lowMemory -- (logical) If we should use the low memory DDA
      %     implementation.  Default: ``false``.
      %
      % Additional parameters passed to class constructor.

      p = inputParser;
      p.addParameter('spacing', 1/20, @isnumeric);
      p.addParameter('polarizability', ...
          @(~, s, r) ott.tmatrix.dda.polarizability.LDR(s, r));
      p.addParameter('relative_index', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct a DDA representation
      dda = ott.tmatrix.dda.Dda.FromShape(...
          'spacing', p.Results.

      % Select the high memory implementation
      if ~p.Results.lowMemory
        dda = ott.tmatrix.dda.DdaHighMem(dda);
      end

      % Calculate the T-matrix
      [varargout{1:nargout}] = ott.tmatrix.Dda(dda, unmatched{:});
    end
  end

  methods
    function [tmatrix, incData, pmData] = Dda(dda_calc, varargin)
      % Construct a T-matrix using a DDA simulation for field calculation
      %
      % Usage
      %   tmatrix = Dda(dda_calc, ...)
      %
      % Parameters
      %   - dda_calc -- (ott.tmatrix.dda.Dda instance) The DDA instance
      %     which will be used for field calculations.
      %
      % Optional named arguments
      %   - Nmax -- (numeric) Size of the VSWF expansion used for the
      %     T-matrix point matching (determines T-matrix rows).
      %     Default: ``ott.utis.ka2nmax(2*pi*shape.maxRadius)`` (may
      %     need different values to give convergence for some shapes).
      %
      %   - ci -- (N numeric) Number of modes to calculate scattering for.
      %     This determines number of T-matrix columns.
      %     Default: ``1:ott.utils.combined_index(Nmax, Nmax)`` (all modes).
      %
      %   - pmradius --(numeric) Radius for point matching.
      %     Default: ``Inf`` (i.e., far-field point matching).
      %
      %   - angulargrid ({theta, phi}) -- Angular grid of points for
      %     calculation of radii.  Default is equally spaced angles with
      %     the number of points determined by Nmax.

      p = inputParser;
      p.addParameter('Nmax', []);
      p.addParameter('pmradius', Inf, ...
          @(x) isnumeric(x) && isscalar(x) && x > 0);
      p.addParameter('ci', []);
      p.addParameter('angulargrid', []);
      p.addParameter('incData', ott.utils.VswfData());
      p.addParameter('pmData', ott.utils.VswfData());
      p.parse(varargin{:});

      % Get or calculate Nmax
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        Nmax = ott.utils.ka2nmax(2*pi*shape.maxRadius);
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
        if dda.zRotSymmetry > 1
          nphi = ceil(nphi ./ dda.zRotSymmetry);
        elseif dda.zRotSymmetry == 0

          % TODO: When we work out how to infinite rotational DDA
          %   we should change this value to 1, for now its 3.
          nphi = 3;
        end

        % Similarly, we can reduce the number of point in theta
        % when we have z-mirror symmetry since we match twice
        if dda.xySymmetry
          ntheta = ceil(ntheta ./ 2);
        end

        theta = ((1:ntheta)-0.5)/ntheta * pi;
        phi = ((1:nphi)-1)/nphi * 2*pi;

        [theta, phi] = meshgrid(theta, phi);

        % We also need to rescale our points by a similar amount
        % This avoids making the problem rank deficient
        if dda.zRotSymmetry > 1
          phi = phi ./ dda.zRotSymmetry;
        end
        if dda.xySymmetry
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

      incData = p.Results.incData;
      pmData = p.Results.pmData;

      % Generate basis set of beams
      Z = sparse(max(ci), numel(ci));
      I = sparse(ci, 1:numel(ci), ones(size(ci)), max(ci), numel(ci));
      vswfBasis = ott.bsc.Bsc([I, Z], [Z, I]);

      % Calculate fields for point matching
      if isfinite(p.Results.pmradius)
        [pmE, pmData] = vswfBasis.efieldRtp(rtp, 'data', pmData);
        pmE = reshape(pmE.vxyz(2:3, :), 3*size(rtp, 2), 2*numel(ci));
      else
        [pmE, pmData] = vswfBasis.efarfield(rtp, 'data', pmData, ...
            'basis', 'outgoing');
        pmE = reshape(pmE.vrtp(2:3, :), 2*size(rtp, 2), 2*numel(ci));
      end

      % Calculate incident fields at dipole locations
      [incE, incData] = vswfBasis.efield(dda.voxels, 'data', incData);
      incE = reshape(incE.vxyz, 3*size(dda.voxels, 2), 2*numel(ci));

      % Get mode indices from ci
      [nmodes, mmodes] = ott.utils.combined_index(ci);

      % Calculate total number of orders
      total_orders_cols = ott.utils.combined_index(max(nmodes), max(nmodes));
      total_orders_rows = ott.utils.combined_index(Nmax, Nmax);

      % Allocate memory for T-matrix data
      data = zeros(2*total_orders_rows, 2*total_orders_cols);

      % rorder is the most expensive dda step, do it first
      for mm = mod(reshape(mmodes, 1, []), 2*dda.zRotSymmetry)

        % TODO: Display progress ((mm+1)/2*z_rotation)

        % Find our indices
        ourIdx = mod(mmodes, 2*dda.zRotSymmetry) == mm;
        ourIdx2 = [ourIdx, ourIdx];

        % Calculate even part
        if z_mirror

          p = mod(mmodes(mm) + nmodes(mm), 2) == 0;
          ourIdxEvn = [p, ~p] & ourIdx2;
          ourIdxOdd = [~p, p] & ourIdx2;

          dipolesEvn = dda.solve(incE(:, ourIdxEvn), ...
              'parity', 'even', 'rorder', mm);
          dipolesOdd = dda.solve(incE(:, ourIdxOdd), ...
              'parity', 'odd', 'rorder', mm);

          % Calculate fields for point matching
          if isfinite(p.Results.pmradius)
            pmEevn = dipolesEvn.efieldRtp(rtp);
            pmEevn = reshape(pmEevn.vxyz, 3*size(rtp, 2), sum(ourIdxEvn));

            pmEodd = dipolesOdd.efieldRtp(rtp);
            pmEodd = reshape(pmEodd.vxyz, 3*size(rtp, 2), sum(ourIdxOdd));
          else
            pmEevn = dipoleEvn.efarfield(rtp);
            pmEevn = reshape(pmEevn.vrtp(2:3, :), ...
                2*size(rtp, 2), sum(ourIdxEvn));

            pmEodd = dipoleOdd.efarfield(rtp);
            pmEodd = reshape(pmEodd.vrtp(2:3, :), ...
                2*size(rtp, 2), sum(ourIdxOdd));
          end

          error('Not yet implemented');

        else
          dipoles = dda.solve(incE(:, ourIdx2), 'rorder', mm);

          % Calculate fields for point matching
          if isfinite(p.Results.pmradius)
            pmEout = dipoles.efieldRtp(rtp);
            pmEout = reshape(pmEout.vxyz, 3*size(rtp, 2), sum(ourIdx2));
          else
            pmEout = dipole.efarfield(rtp);
            pmEout = reshape(pmEout.vrtp(2:3, :), 2*size(rtp, 2), sum(ourIdx2));
          end

          if z_rotation == 1
            % No rotational symmetry
            data(:, [ci, total_orders_cols+ci]) = pmE \ pmEout;
          else

            [alln, allm] = ott.utils.combined_index((1:total_orders_rows).');

            if z_rotation > 1
              % Calculate which modes preseve symmetry, m = +/- ip
              axial_modes = mod(allm - mmodes(mm), z_rotation) == 0;
            elseif z_rotation == 0
              % Modes only scatter to modes with same m
              axial_modes = allm == mmodes(mm);
            else
              error('Invalid z_rotation value');
            end

            modes = [axial_modes; axial_modes];

            pq1 = pmE(:, modes) \ pmEout(:, [ourIdx, 0*ourIdx]);
            pq2 = pmE(:, modes) \ pmEout(:, [0*ourIdx, ourIdx]);

            data(modes, ci) = pq1;
            data(modes, total_orders_cols+ci) = pq2;
          end
        end
      end

      % Package T-matrix
      tmatrix.data = data;
      tmatrix = tmatrix.setType('scattered');
    end
  end
end
