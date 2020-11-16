classdef LgParaxial < ott.bsc.Bsc
% Implements point matching using a reduced LG basis in the paraxial limit.
% Inherits from :class:`Bsc`.
% Infernally uses :meth:`Pointmatch.FromFarfield`.
%
% This class uses an LG basis for point matching.  Instead of using a
% 2-D grid of points over the far-field each LG radial mode can be
% matched with a 1-D array of points, greatly reducing the number of
% points needed and consequently reducing the computation time.
%
% Properties
%   - a         -- Beam shape coefficients `a` vector
%   - b         -- Beam shape coefficients `b` vector
%   - waist     -- Paraxial beam waist
%   - lmode     -- Azimuthal mode indices
%   - pmode     -- Radial mode indices (>= 0)
%   - polmode   -- Polarisation mode indices (+/- 1)
%   - mapping   -- Paraxial coordinate mapping (enum)
%   - truncation_angle -- Far-field truncation angle [0, pi]
%   - Nmax      -- (Dependent) Truncation number for VSWF coefficients
%   - power     -- (Dependent) Power of the beam shape coefficients
%
% Static methods
%   - FromIgMode     -- Construct beam from Ince-Gaussian mode
%   - FromLgMode     -- Construct beam from Laguerre-Gaussian mode
%   - FromHgMode     -- Construct beam from Hermite-Gaussian mode
%   - WaistFromNa    -- Convert numerical aperture to paraxial waist
%   - ParaxialWaist  -- Compute paraxial beam waist for high order mode

% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    waist      % Paraxial beam waist
    lmode      % Azimuthal mode indices
    pmode      % Radial mode indices (>= 0)
    polmode    % Polarisation mode indices (+/- 1)
    mapping    % Paraxial coordinate mapping (enum)
    truncation_angle  % Far-field truncation angle [0, pi]
  end

  methods (Static)
    function [bsc, weights, data] = FromIgMode(waist, lmode, porder, parity, ...
          ellipticity, polbasis, polfield, varargin)
      % Construct a Ince-Gaussian paraxial mode
      %
      % Usage
      %   bsc = LgParaxialBasis.FromIgMode(waist, lmode, porder, parity,
      %   ellipticity, polbasis, polfield)
      %
      %   [bsc, weights] = FromIgMode(...) -- As above, but returns the
      %   LgParaxial beam and a vector of weights.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [m].
      %
      %   - lmode (numeric) -- Azimuthal mode order.
      %
      %   - porder (numeric) -- Paraxial mode order.
      %
      %   - ellipticity (numeric) -- Ellipticity of the beam.
      %
      %   - parity (enum) -- Parity.  Either 'even' or 'odd'.
      %
      %   - polbasis (enum) -- Polarisation basis.  Either 'polar' or
      %     'cartesian'.
      %
      %   - polfield (2 numeric) -- Field in either the theta/phi or
      %     x/y directions (depending on basis).
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
      assert(sum(strcmpi(polbasis, {'cartesian', 'polar'})) == 1, ...
          'polbasis must be ''cartesian'' or ''polar''');
      assert(isnumeric(polfield) && numel(polfield) == 2, ...
          'polfield should be 2 element numeric vector');

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

      weights = [1, -1i; 1, 1i] * polfield(:);
      weights = repelem(weights, nummodes, 1);
      weights = repmat(modeweights(row, keepz).', 2, 1) .* weights;

      % Generate beam
      [bsc, data] = ott.bsc.LgParaxial(waist, lmode, pmode, polmode, ...
          varargin{:});

      if nargout == 1
        bsc = bsc * weights;
      end
    end

    function [bsc, weights, vswfData] = FromLgMode(waist, lmode, ...
        pmode, polbasis, polfield, varargin)
      % Construct a Laguerre-Gaussian paraxial mode
      %
      % Usage
      %   bsc = LgParaxialBasis.FromLgMode(waist, lmode, pmode,
      %   polbasis, polfield)
      %
      %   [bsc, weights] = FromLgMode(...) -- As above, but returns the
      %   LgParaxial beam and a vector of weights.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [m].
      %
      %   - lmode (numeric) -- Azimuthal mode number.
      %
      %   - pmode (numeric) -- Radial mode number.
      %
      %   - polbasis (enum) -- Polarisation basis.  Either 'polar' or
      %     'cartesian'.
      %
      %   - polfield (2 numeric) -- Field in either the theta/phi or
      %     x/y directions (depending on basis).
      %
      % Unmatched arguments are passed to constructor.

      % Check parameters
      assert(isnumeric(waist) && isscalar(waist), ...
          'waist must be numeric scalar');
      assert(isnumeric(lmode) && isscalar(lmode), ...
          'lmode must be numeric scalar');
      assert(isnumeric(pmode) && isscalar(pmode), ...
          'pmode must be numeric scalar');
      assert(sum(strcmpi(polbasis, {'cartesian', 'polar'})) == 1, ...
          'polbasis must be ''cartesian'' or ''polar''');
      assert(isnumeric(polfield) && numel(polfield) == 2, ...
          'polfield should be 2 element numeric vector');

      lmode = [lmode, lmode];
      pmode = [pmode, pmode];
      polmode = [-1, 1];
      weights = [1, -1i; 1, 1i] * polfield(:);

      % Generate beam
      [bsc, vswfData] = ott.bsc.LgParaxial(waist, lmode, pmode, polmode, ...
          varargin{:});

      if nargout == 1
        bsc = bsc * weights;
      end
    end

    function [bsc, weights, vswfData] = FromHgMode(waist, mmode, nmode, ...
        polbasis, polfield, varargin)
      % Construct a Hermite-Gaussian paraxial mode
      %
      % Usage
      %   bsc = LgParaxialBasis.FromHgMode(waist, mmode, nmode,
      %       polbasis, polfield, ...)
      %
      %   [bsc, weights, data] = FromHgMode(...) -- As above, but returns the
      %   LgParaxial beam, a vector of weights, and VSWF data for re-use.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [m].
      %
      %   - mmode (numeric) -- Mode number.
      %
      %   - nmode (numeric) -- Mode number.
      %
      %   - polbasis (enum) -- Polarisation basis.  Either 'polar' or
      %     'cartesian'.
      %
      %   - polfield (2 numeric) -- Field in either the theta/phi or
      %     x/y directions (depending on basis).
      %
      % Unmatched arguments are passed to constructor.

      % Check parameters
      assert(isnumeric(waist) && isscalar(waist), ...
          'waist must be numeric scalar');
      assert(isnumeric(mmode) && isscalar(mmode) && mmode >= 0, ...
          'mmode must be positive numeric scalar');
      assert(isnumeric(nmode) && isscalar(nmode) && nmode >= 0, ...
          'nmode must be positive numeric scalar');
      assert(sum(strcmpi(polbasis, {'cartesian', 'polar'})) == 1, ...
          'polbasis must be ''cartesian'' or ''polar''');
      assert(isnumeric(polfield) && numel(polfield) == 2, ...
          'polfield should be 2 element numeric vector');

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

      weights = [1, -1i; 1, 1i] * polfield(:);
      weights = repelem(weights, nummodes, 1);
      weights = repmat(modeweights(row, keepz).', 2, 1) .* weights;

      % Generate beam
      [bsc, vswfData] = ott.bsc.LgParaxial(waist, lmode, pmode, polmode, ...
          varargin{:});

      if nargout == 1
        bsc = bsc * weights;
      end
    end

    function w = ParaxialWaist(paraxial_order)
      % Computes the paraxial beam waist for a high-order Gaussian beam.
      %
      % Based on the OTTv1 function ``paraxial_beam_waist``.
      %
      % Usage
      %   w0 = ParaxialWaist(paraxial_order)

      w = 1.; %Beam waist in normalized units.

      if (paraxial_order ~= 0)
        invL=1./abs(paraxial_order );
        zz = exp(-(abs(paraxial_order )+2.)*invL);

        % Choose a good starting guess (this converges within 3 iterations
        % for l between 1 and 10000+
        w=-(1.+2*sqrt(invL)+invL);

        w0=-w;

        while (abs(w-w0)>0.00001)
          w0=w;
          expw = exp(w);

          % Use Newton's method to find the root
          % Would typically find the real root, ...
          % This finds the W_{-1}(z) local to the beam waist of an LG beam.
          w=w0-(w0*expw+zz)/(expw+w0*expw);

        end

        %Beam waist in normalized units
        w = sqrt(-abs(paraxial_order )/2.*w); 

      end
    end

    function waist = WaistFromNa(na, index_medium, paraxial_order, mapping)
      % Convert numerical aperture to paraxial waist
      %
      % Usage
      %   waist = ott.bsc.LgParaxial.WaistFromNa(NA, ...
      %   index_medium, porder, mapping)
      %
      % Parameters
      %   - NA (numeric) -- Numerical aperture
      %
      %   - index_medium (numeric) -- Refractive index in medium.
      %
      %   - porder (numeric) -- Paraxial mode order.
      %
      %   - mapping (enum) -- Paraxial mapping for back aperture.
      %     Can be either 'tan' or 'sin'.

      % Get beam angle
      angle_rad = ott.utils.na2angle(na, index_medium);

      % Calculate paraxial waist
      w0 = ott.bsc.LgParaxial.ParaxialWaist(paraxial_order);

      % Calculate scaling factor
      switch mapping
        case 'sin'
          wscaling = 1/sin(abs(angle_rad));
        case 'tan'
          wscaling = 1/tan(abs(angle_rad));
        otherwise
          error('Unknown mapping option');
      end

      waist = w0 .* wscaling;
    end
	end

  methods
    function [bsc, vswfData] = LgParaxial(varargin)
      % Generate a LG paraxial mode
      %
      % Usage
      %   [bsc, vswfData] = LgParaxial(waist, lmode, pmode, polmod, ...)
      %
      % Parameters
      %   - waist (numeric) -- Paraxial beam waist
      %   - lmode (numeric) -- Azimuthal mode indices
      %   - pmode (numeric) -- Radial mode indices
      %   - polmode (numeric) -- Polarisation mode (+/- 1)
      %
      % Named parameters
      %   - mapping (enum) -- Mapping method for paraxial far-field.
      %     Can be either 'sin' or 'tan' (small angle).
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).  Default: ``'sin'``.
      %
      %   - truncation_angle (numeric) -- Truncation angle (in radians)
      %     for input aperture cut-off.  Default: ``pi/2``,
      %     corresponding to one full hemisphere.
      %
      %   - Nmax (numeric) -- Nmax for calculation.  Default: ``100``
      %     (provides good approximation for most tightly focussed beams).
      %
      %   - data (ott.utils.VswfData) -- Data for repeated beam calculation.
      %     Default: ``ott.utils.VswfData()``.

      p = inputParser;
      p.addOptional('waist', [], @isnumeric);
      p.addOptional('lmode', [], @isnumeric);
      p.addOptional('pmode', [], @isnumeric);
      p.addOptional('polmode', [], @isnumeric);
      p.addParameter('mapping', 'sin');
      p.addParameter('truncation_angle', pi/2);
      p.addParameter('data', ott.utils.VswfData());
      p.addParameter('Nmax', 100);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bsc = bsc@ott.bsc.Bsc(unmatched{:});
      bsc.waist = p.Results.waist;
      bsc.mapping = p.Results.mapping;
      bsc.truncation_angle = p.Results.truncation_angle;

      lmode = p.Results.lmode;
      pmode = p.Results.pmode;
      polmode = p.Results.polmode;
      
      % Check waist range
      if bsc.waist < 1e-2
        warning('ott:bsc:LgParaxial:small_waist', ...
          'Requested waist is less than 1%% of wavelength');
      end

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

      % Get Nmax (if not already set)
      Nmax = p.Results.Nmax;

       % Get the mode indices for this Nmax
      total_modes = ott.utils.combined_index(Nmax, Nmax);
      ci = 1:total_modes;
      [~,mm] = ott.utils.combined_index(ci');

      % Calculate theta points and angular grid
      ntheta = 2*(Nmax + 1);
      nphi = 1;
      [theta, phi] = ott.utils.angulargrid(ntheta, nphi);
      rtp = [theta(:), phi(:)].';

      newmodes = ott.bsc.Bsc.empty();
      sorder = [];
      
      vswfData = p.Results.data;

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
        Ertp = bsc.paraxial_fields(theta, azimuthal_modes, ...
            radial_modes, polarisation_modes);

        % Construct a beam with far-field point matching
        [modes, vswfData] = ott.bsc.Bsc.PmFarfield(rtp, Ertp, ci(midx), ...
            'data', vswfData);
        newmodes = [newmodes, modes.split()];
      end

      % Sort new rows (to match order of inputs)
      newmodes(1, sorder) = newmodes(1, :);

      % Allocate results
      bsc = repmat(bsc, 1, Nmodes);
      bsc = bsc.setCoefficients(newmodes);
      for ii = 1:numel(bsc)
        bsc(ii).lmode = lmode(ii);
        bsc(ii).pmode = pmode(ii);
        bsc(ii).polmode = polmode(ii);
      end
    end
  end
  
  methods (Hidden)
    function [rw, dr] = paraxial_scaling(bsc, theta)
      % Calculate paraxial scaling factors

      owaist = bsc.waist;

      switch bsc.mapping
        case 'sin'
          rw = 2*(owaist).^2 .* sin(theta).^2;
          dr = owaist .* abs(cos(theta));
        case 'tan'
          rw = 2*(owaist).^2 .* tan(theta).^2;
          dr = owaist .* (sec(theta)).^2;
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

      Efield = zeros(2, numel(theta), numel(azimuthal_mode));

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
        Efield(1, :, ii) = Etheta;
        Efield(2, :, ii) = Ephi;

      end

      % Apply truncation angle
      mask = theta < pi - bsc.truncation_angle;
      Efield(:, mask, :) = 0.0;
    end
  end
  
  methods % Getters/setters

    function bsc = set.mapping(bsc, val)
      assert(any(strcmpi(val, {'sin', 'tan'})), ...
          'mapping must be ''sin'' or ''tan''');
      bsc.mapping = val;
    end

    function bsc = set.truncation_angle(bsc, val)
      assert(isnumeric(val) && isscalar(val), ...
        'truncation_angle must be numeric scalar');
      assert(val >= 0.0 && val <= pi, ...
          'truncation_angle must be between 0 and pi');
      bsc.truncation_angle = val;
    end

    function bsc = set.waist(bsc, val)
      assert(isnumeric(val) && isscalar(val), ...
        'waist must be numeric scalar');
      bsc.waist = val;
    end

    function bsc = set.lmode(bsc, val)
      assert(isnumeric(val) && isscalar(val), ...
        'lmode must be numeric scalar');
      bsc.lmode = val;
    end

    function bsc = set.pmode(bsc, val)
      assert(isnumeric(val) && isscalar(val), ...
        'pmode must be numeric scalar');
      bsc.pmode = val;
    end

    function bsc = set.polmode(bsc, val)
      assert(isnumeric(val) && isscalar(val), ...
        'polmode must be numeric scalar');
      bsc.polmode = val;
    end
  end
end

