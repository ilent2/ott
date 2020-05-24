classdef LgParaxialBasis < ott.beam.vswf.Bsc
% Implements point matching using a reduced LG basis in the paraxial limit.
% Inherits from :class:`ott.beam.vswf.Bsc`, uses
% :class:`ott.beam.vswf.FarfieldPm` for point-matching.
%
% This class uses an LG basis for point matching.  Instead of using a
% 2-D grid of points over the far-field each LG radial mode can be
% matched with a 1-D array of points, greatly reducing the number of
% points needed and consequently reducing the computation time.
%
% Properties
%   - a, b          -- Beam shape coefficients.
%   - lmode         -- Azimuthal mode indices
%   - pmode         -- Radial mode indices
%   - polmode       -- Polarisation mode indices (+/-1)
%   - mapping           -- Paraxial coordinate mapping
%   - truncation_angle  -- Far-field truncation angle [0, pi]
%
% Methods
%   - addModes       -- Add LG modes to the beam array
%
% Static methods
%   - FromIgMode     -- Construct beam from Ince-Gaussian mode
%   - FromLgMode     -- Construct beam from Laguerre-Gaussian mode
%   - FromHgMode     -- Construct beam from Hermite-Gaussian mode
%   - FromParaxial   -- Construct beam from LG decomposition of paraxial field.
%   - paraxial_weights -- Calculates the weights vector for a paraxial field.
%   - empty            -- Construct an empty beam array.
%
% Hidden methods
%   - paraxial_fields   -- Used by `addModes`
%   - paraxial_scaling  -- Used by paraxial_fields

% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    waist             % Paraxial beam waist
    lmode             % Azimuthal mode indices
    pmode             % Radial mode indices
    polmode           % Polarisation mode (+/- 1)
    mapping           % Paraxial coordinate mapping
    truncation_angle  % Far-field truncation angle [0, pi]
  end

  methods (Static)
    function [lmodes, pmodes, polmodes, weights] = paraxial_weights(...
          waist, Pmax, xy, Exy, varargin)
      % Performs point-matching to calculate LG mode weights
      %
      % Performs point-matching twice, first for the left circular
      % polarisation basis and then for the right basis.
      %
      % Usage
      %   [lmode, pmode, polmode, weights] = paraxial_weights(
      %       waist, Pmax, xy, Exy, ...)
      %
      % Parameters
      %   - waist (numeric) -- Scaling factor for paraxial waist.
      %   - Pmax (numeric) -- Maximum paraxial order.  This corresponds
      %     to :math:`2*p_{mode} + \abs(l_{mode}) \leq P_{max}`.
      %   - xy (2xN numeric) -- Paraxial coordinates of samples.
      %   - Exy (2xN numeric) -- Sampled field to point-match.
      %
      % Returns
      %   - lmodes (numeric) -- Azimuthal mode indices.
      %   - pmodes (numeric) -- Radial mode indices.
      %   - weights (numeric) -- Weights of matched modes.
      %
      % Optional named arguments
      %   - filter_tol (numeric) -- Minimum power to keep modes.
      %     Default: ``1.0e-6*max(abs(Exy))``.

      p = inputParser;
      p.addParameter('filter_tol', 1.0e-6*max(abs(Exy)));
      p.parse(varargin{:});

      assert(isnumeric(xy) && size(xy, 1) == 2, 'xy must be 2xN numeric');
      assert(isnumeric(Exy) && size(Exy, 1) == 2, 'Exy must be 2xN numeric');
      assert(size(xy, 2) == size(Exy, 2), 'xy and Exy must have same length');
      assert(isnumeric(Pmax) && isscalar(Pmax), ...
          'Pmax must be numeric scalar');
      assert(isnumeric(waist) && isscalar(waist), ...
          'waist must be numeric scalar');

      % Convert from xy polarisation to lr
      Elr = [1, -1i; 1, 1i] * Exy;

      % Get polar coordinates
      rho = vecnorm(xy) ./ waist;
      phi = atan2(xy(2, :), xy(1, :));

      radial_modes = 0:Pmax/2;
      azimuthal_max = floor(Pmax - 2*radial_modes);
      Nmodes = sum(2*azimuthal_max + 1);

      % Generate coefficient matrix
      coefficient_matrix = zeros(size(xy, 2), Nmodes);
      lmodes = zeros(Nmodes, 1);
      pmodes = zeros(Nmodes, 1);
      idx = 1;

      for ii = 1:length(radial_modes)

        rmode = radial_modes(ii);
        azimuthal_modes = -azimuthal_max(ii):azimuthal_max(ii);

        % Store mode indices
        lmodes(idx:idx+numel(azimuthal_modes)) = azimuthal_modes;
        pmodes(idx:idx+numel(azimuthal_modes)) = rmode;

        for jj = 1:length(azimuthal_modes)

          amode = azimuthal_modes(jj);

          % Generate LG modes
          L = ot.utils.laguerre(rmode,abs(amode), rho) .* exp(1i*phi);
          coefficient_matrix(:, idx+jj-1) = L(:);
        end

        % Increment array index
        idx = idx + numel(azimuthal_modes);
      end

      % Do point matching
      weights = coefficient_matrix \ Elr.';

      % Find non-zero weights
      keep = abs(weights) >= p.Results.keep_top;

      % Package outputs
      polmodes = [-1*ones(sum(keep(:, 1)), 1); ones(sum(keep(:, 2)), 1)];
      lmodes = [lmodes(keep(:, 1)); lmodes(keep(:, 2))];
      pmodes = [pmodes(keep(:, 1)); pmodes(keep(:, 2))];
      weights = weights(keep);

    end

    function bsc = FromParaxial(waist, Pmax, xy, Exy, varargin)
      % Construct a beam from LG decomposition of paraxial field.
      %
      % Only calculates bsc-modes for LG-modes with significant power.
      %
      % Optimal LG-mode basis size depends on the overlap between the
      % target field and the LG basis scaling factor (`waist`).
      % For small LG-mode bases, this method may be faster than using
      % a plane wave basis; for arbitrary fields, it is probably better
      % to use a plane wave or Bessel basis instead.
      %
      % Usage
      %   bsc = LgParaxialBasis.FromParaxial(waist, Pmax, xy, Exy, ...)
      %
      % Parameters
      %   - waist (numeric) -- Waist for the LG-mode basis.
      %   - Pmax (numeric) -- Maximum LG-mode paraxial order.
      %   - xy (2xM numeric) -- Coordinates for point matching.
      %   - Exy (2xM numeric) -- Paraxial far-field for point matching.
      %
      % Optional named arguments
      %   - filter_tol (numeric) -- Minimum power to keep modes.
      %     Default: ``1.0e-6*max(abs(Exy))``.
      %
      % Unmatched arguments are passed to class constructor.

      p = inputParser;
      p.addParameter('filter_tol', 1.0e-6*max(abs(Exy)));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Calculate weights
      [lmode, pmode, polmode, weights] = ...
          ott.beam.vswf.LgParaxialBasis.paraxial_weights(...
          waist, Pmax, xy, Exy, 'filter_tol', p.Results.filter_tol);

      % Generate beam
      bsc = ott.beam.vswf.LgParaxialBasis(waist, lmode, pmode, polmode, ...
          'weights', weights, unmatched{:});
    end

    function bsc = FromIgMode(waist, lmode, porder, parity, ...
          ellipticity, polarisation, varargin)
      % Construct a Ince-Gaussian paraxial mode
      %
      % Usage
      %   bsc = LgParaxialBasis.FromIgMode(waist, lmode, porder, parity,
      %   ellipticity, polarisation)
      %
      % Unmatched arguments are passed to constructor.

      % Check parameters
      assert(isnumeric(waist) && isscalar(waist), ...
          'waist must be numeric scalar');
      assert(isnumeric(lmode) && isscalar(lmode), ...
          'lmode must be numeric scalar');
      assert(isnumeric(porder) && isscalar(porder), ...
          'pmode must be numeric scalar');
      assert(sum(strcmpi(parity, {'even', 'odd'})) == 1, ...
          'parity must be ''even'' or ''odd''');
      assert(isnumeric(ellipticity) && isscalar(ellipticity), ...
          'ellipticity must be numeric scalar');
      assert(isnumeric(polarisation) && numel(polarisation) == 2, ...
          'polarisation should be 2 element numeric vector');

      [modeweights,initial_mode,final_mode] = ...
         ott.utils.paraxial_transformation_matrix(porder,0,[2,ellipticity],0);

      if strcmpi(parity, 'odd')
        parity = 0;
      else
        parity = 1;
      end

      [row]=find(and(final_mode(:,2)==lmode, ...
          final_mode(:,3)==parity),1);

      if porder > 1 && isempty(row)
        error('invalid parameters to satisfy parity conventions');
      end

      keepz = (abs(modeweights(row, :)) > 0);
      nummodes = sum(keepz);

      lmode = initial_mode(keepz, 2);
      pmode = initial_mode(keepz, 1);

      polmode = repelem([-1; 1], nummodes, 1);
      lmode = [lmode; lmode];
      pmode = [pmode; pmode];

      weights = [1, -1i; 1, 1i] * polarisation(:);
      weights = repelem(weights, nummodes, 1);
      weights = repmat(modeweights(row, keepz).', 2, 1) .* weights;

      % Generate beam
      bsc = ott.beam.vswf.LgParaxialBasis(waist, lmode, pmode, polmode, ...
          'weights', weights, varargin{:});
    end

    function bsc = FromLgMode(waist, lmode, pmode, polarisation, varargin)
      % Construct a Laguerre-Gaussian paraxial mode
      %
      % Usage
      %   bsc = LgParaxialBasis.FromLgMode(waist, lmode, pmode, polarisation)
      %
      % Parameters
      %   - lmode (numeric) -- Azimuthal mode
      %   - pmode (numeric) -- Radial mode
      %   - polarisation (2-numeric) -- Jones vector for polarisation.
      %
      % Unmatched arguments are passed to constructor.

      % Check parameters
      assert(isnumeric(waist) && isscalar(waist), ...
          'waist must be numeric scalar');
      assert(isnumeric(lmode) && isscalar(lmode), ...
          'lmode must be numeric scalar');
      assert(isnumeric(pmode) && isscalar(pmode), ...
          'pmode must be numeric scalar');
      assert(isnumeric(polarisation) && numel(polarisation) == 2, ...
          'polarisation should be 2 element numeric vector');

      lmode = [lmode, lmode];
      pmode = [pmode, pmode];
      polmode = [-1, 1];
      weights = [1, -1i; 1, 1i] * polarisation(:);

      % Generate beam
      bsc = ott.beam.vswf.LgParaxialBasis(waist, lmode, pmode, polmode, ...
          'weights', weights, varargin{:});
    end

    function bsc = FromHgMode(waist, mmode, nmode, polarisation, varargin)
      % Construct a Hermite-Gaussian paraxial mode
      %
      % Usage
      %   bsc = LgParaxialBasis.FromHgMode(waist, mmode, nmode,
      %       polarisation, ...)
      %
      % Unmatched arguments are passed to constructor.

      % Check parameters
      assert(isnumeric(waist) && isscalar(waist), ...
          'waist must be numeric scalar');
      assert(isnumeric(mmode) && isscalar(mmode), ...
          'mmode must be numeric scalar');
      assert(isnumeric(nmode) && isscalar(nmode), ...
          'nmode must be numeric scalar');
      assert(isnumeric(polarisation) && numel(polarisation) == 2, ...
          'polarisation should be 2 element numeric vector');

      paraxial_order=nmode+mmode;
      [modeweights,initial_mode,final_mode] = ...
          ott.utils.paraxial_transformation_matrix(paraxial_order,0,1,0);
      [row] = find(final_mode(:,1)==mmode,1);

      keepz = (abs(modeweights(row, :)) > 0);
      nummodes = sum(keepz);

      lmode = initial_mode(keepz, 2);
      pmode = initial_mode(keepz, 1);

      polmode = repelem([-1; 1], nummodes, 1);
      lmode = [lmode; lmode];
      pmode = [pmode; pmode];

      weights = [1, -1i; 1, 1i] * polarisation(:);
      weights = repelem(weights, nummodes, 1);
      weights = repmat(modeweights(row, keepz).', 2, 1) .* weights;

      % Generate beam
      bsc = ott.beam.vswf.LgParaxialBasis(waist, lmode, pmode, polmode, ...
          'weights', weights, varargin{:});
    end

    function bsc = empty(varargin)
      % Create an empty LgParaxialBasis beam collection
      %
      % Usage
      %   beam = LgParaxialBasis.empty()
      %
      % Additional arguments are passed to constructor.

      bsc = ott.beam.vswf.LgParaxialBasis(1.0, [], [], [], varargin{:});
    end
  end

  methods
    function bsc = LgParaxialBasis(varargin)
      % Generate a LG mode basis
      %
      % Usage
      %   bsc = LgParaxialBasis(waist, lmode, pmode, polmod, ...)
      %
      %   bsc = LgParaxialBasis(waist, ...)
      %   Constructs a beam array with not beams.  Modes can be
      %   added after construction using :meth:`addModes`.
      %
      %   bsc = LgParaxialBasis(...)
      %   Constructs a empty beam array with default parameters.
      %
      % Parameters
      %   - waist (numeric) -- Paraxial beam waist.  Default: ``1.0``.
      %   - lmode (numeric) -- Azimuthal mode indices
      %   - pmode (numeric) -- Radial mode indices
      %   - polmode (numeric) -- Polarisation mode (+/- 1)
      %
      % Named parameters
      %   - weights (numeric) -- Weights to apply to each beam.
      %     Default: ``[]`` (i.e., don't apply weights).  Weights can
      %     also be applied by multiplying by a weights matrix,
      %     for example ``single_bsc = bsc * weights`` or
      %     ``multiple_bsc = bsc .* weights``.
      %
      %   - mapping (enum) -- Mapping method for paraxial far-field.
      %     Can be either 'sintheta' or 'tantheta' (small angle).
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).  Default: ``'sintheta'``.
      %
      %   - truncation_angle (numeric) -- Truncation angle (in radians)
      %     for input aperture cut-off.  Default: ``pi/2``,
      %     corresponding to one full hemisphere.

      p = inputParser;
      p.addOptional('arg1', [], @isnumeric);
      p.addOptional('arg2', [], @isnumeric);
      p.addOptional('arg3', [], @isnumeric);
      p.addOptional('arg4', [], @isnumeric);
      p.addParameter('weights', []);
      p.addParameter('mapping', 'sintheta');
      p.addParameter('truncation_angle', pi/2);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      % Construct base
      unmatched = ott.utils.unmatchedArgs(p);
      bsc = bsc@ott.beam.vswf.Bsc(unmatched{:});

      arg1 = p.Results.arg1;
      arg2 = p.Results.arg2;
      arg3 = p.Results.arg3;
      arg4 = p.Results.arg4;

      % Check number of arguments
      num_args = ~isempty(arg1) + ~isempty(arg2) + ~isempty(arg3) + ~isempty(arg4);
      assert(num_args == 0 || num_args == 1 || num_args == 4, ...
          'Must provide either 0, 1 or 4 positional arguments');

      % Store paraxial waist
      if num_args >= 1
        bsc.waist = arg1;
      else
        bsc.waist = 1.0;
      end

      % Store paraxial mapping parameter and truncation_angle
      bsc.mapping = p.Results.mapping;
      bsc.truncation_angle = p.Results.truncation_angle;

      % Add modes
      if num_args == 4
        bsc = bsc.addModes(arg2, arg3, arg4);
      end

      % Apply weights
      if ~isempty(p.Results.weights)
        bsc = bsc * p.Results.weights;
      end
    end

    function bsc = addModes(bsc, lmode, pmode, polmode)
      % Add modes to the beam array
      %
      % Usage
      %   bsc = bsc.addModes(lmode, pmode, polmode)
      %
      % See constructor for parameter information.

      Nmodes = max([numel(lmode), numel(pmode), numel(polmode)]);

      % Check sizes
      assert(numel(lmode) == Nmodes || numel(lmode) == 1, ...
        'lmode must be scalar or match size of pmode/polmode');
      assert(numel(pmode) == Nmodes || numel(pmode) == 1, ...
        'pmode must be scalar or match size of lmode/polmode');
      assert(numel(polmode) == Nmodes || numel(polmode) == 1, ...
        'polmode must be scalar or match size of lmode/pmode');

      % Check value ranges
      assert(all(round(lmode) == lmode), ...
        'ott:vswf:LaguerreGaussian:invalid_azimuthal_mode', ...
        'lmode must be integers');
      assert(all(pmode >= 0 & round(pmode) == pmode), ...
        'ott:vswf:LaguerreGaussian:invalid_radial_mode', ...
        'pmode must be positive integers');
      assert(all(polmode == -1 | polmode == 1), ...
          'polmode must be +/- 1');

      % Grow scalar inputs (if required)
      if numel(lmode) == 1, lmode = repmat(lmode, 1, Nmodes); end
      if numel(pmode) == 1, pmode = repmat(pmode, 1, Nmodes); end
      if numel(polmode) == 1, polmode = repmat(polmode, 1, Nmodes); end

      % Find which m-columns we need to solve for
      um = unique(lmode(:) + polmode(:));
      
      % Default Nmax (if not already set)
      if bsc.Nmax == 0
        bsc.Nmax = 100;
      end

      % Get the mode indices for this Nmax
      total_modes = ott.utils.combined_index(bsc.Nmax, bsc.Nmax);
      [nn,mm] = ott.utils.combined_index((1:total_modes)');

      % Calculate theta points and angular grid
      ntheta = 2*(bsc.Nmax + 1);
      nphi = 1;
      [theta, phi] = ott.utils.angulargrid(ntheta, nphi);
      locations = [theta(:), phi(:)].';

      newmodes = ott.beam.vswf.FarfieldPm.empty();
      sorder = [];

      % Do point matching of each umodes together
      for m = um.'

        % Get LG modes for current iteration
        idx = m == (lmode(:) + polmode(:));
        azimuthal_modes = lmode(idx);
        radial_modes = pmode(idx);
        polarisation_modes = polmode(idx);

        % Get BSC modes for current iteration
        midx = m == mm;

        % Store indices (for sorting later)
        sorder = [sorder; find(idx)];

        % Construct field pattern
        Efield = bsc.paraxial_fields(theta, azimuthal_modes, ...
            radial_modes, polarisation_modes);

        % Construct a beam with far-field point matching
        newmodes = [newmodes, ott.beam.vswf.FarfieldPm(...
            nn(midx), mm(midx), locations, Efield)];
      end

      % Sort new rows (to match order of inputs)
      newmodes(1, sorder) = newmodes(1, :);

      % Append bsc data
      bsc = bsc.setCoefficients([bsc.a, newmodes.a], [bsc.b, newmodes.b]);
      bsc.lmode = [bsc.lmode, lmode];
      bsc.pmode = [bsc.pmode, pmode];
      bsc.polmode = [bsc.polmode, polmode];

    end
  end

  methods (Hidden)
    function [rw, dr] = paraxial_scaling(bsc, theta)
      % Calculate paraxial scaling factors

      waist = bsc.waist;

      switch bsc.mapping
        case 'sintheta'
          rw = 2*(waist).^2 .* sin(theta).^2;
          dr = waist .* abs(cos(theta));
        case 'tantheta'
          rw = 2*(waist).^2 .* tan(theta).^2;
          dr = waist .* (sec(theta)).^2;
        otherwise
          error('Unknown mapping option');
      end

    end

    function Efield = paraxial_fields(bsc, theta, ...
        azimuthal_mode, radial_mode, polarisation)
      % Construct paraxial field slices for specified LG mode
      %
      % Base on code by Alexander Stilgoe

      [rw, dr] = bsc.paraxial_scaling(theta);

      Efield = zeros(2*numel(theta), numel(azimuthal_mode));

      for ii = 1:numel(azimuthal_mode)

        norm_paraxial=sqrt(2*factorial(radial_mode(ii)) ...
            /(pi*factorial(radial_mode(ii)+abs(azimuthal_mode(ii)))));
        L = ott.utils.laguerre(radial_mode(ii),abs(azimuthal_mode(ii)),rw);
        beam_envelope = norm_paraxial ...
          .*rw.^abs(azimuthal_mode(ii)/2) .* L ...
          .* exp(-rw/2 + 1i*pi/2 ...
          *(radial_mode(ii)*2+abs(azimuthal_mode(ii))+1));
        mode_input_power=sqrt(sum(2*pi*...
            abs(beam_envelope).^2.*sqrt(rw/2).*abs(dr)));
        aperture_power_normalization= ...
            sqrt(sum(2*pi*abs(beam_envelope).^2.*sin(theta)));

        beam_envelope= beam_envelope ...
            ./aperture_power_normalization*mode_input_power;

        % Apply polarisation (XY Etheta/Ephi <-> Exy when phi = 0)
        Etheta = -1 .* beam_envelope;
        Ephi = 1i .* beam_envelope .* polarisation(ii);

        % Package output
        Efield(:, ii) = [Etheta; Ephi];

      end

      % Apply truncation angle
      mask = theta < pi - bsc.truncation_angle;
      Efield([mask; mask], :) = 0.0;
    end
  end

  methods % Getters/setters
    % mapping  -- set the paraxial mapping function
    % truncation_angle  -- set the truncation angle
    % waist  -- set the beam waist

    function bsc = set.mapping(bsc, val)
      assert(any(strcmpi(val, {'sintheta', 'tantheta'})), ...
          'mapping must be ''sintheta'' or ''tantheta''');
      assert(numel(bsc) == 0, ...
          'Can''t change mapping once beams have been added');
      bsc.mapping = val;
    end

    function bsc = set.truncation_angle(bsc, val)
      assert(isnumeric(val) && isscalar(val), ...
        'truncation_angle must be numeric scalar');
      assert(val >= 0.0 && val <= pi, ...
          'truncation_angle must be between 0 and pi');
      assert(numel(bsc) == 0, ...
          'Can''t change truncation_angle once beams have been added');
      bsc.truncation_angle = val;
    end

    function bsc = set.waist(bsc, val)
      assert(isnumeric(val) && isscalar(val), ...
        'waist must be numeric scalar');
      assert(numel(bsc) == 0, ...
          'Can''t change waist once beams have been added');
      bsc.waist = val;
    end
  end
end
