classdef Annular < ott.bsc.Bsc
% Bsc specialisation for Bessel-like beams
%
% Bessel beams and other annular beams can be translated axially with
% only a phase shift applied to the beam shape coefficients.
% This class overloads the axial translation function to implement this.
%
% Properties
%   - theta       -- Angle describing annular
%
% Static methods
%   - FromBessel  -- Construct Annular beam from Bessel specification.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    theta         % Angle describing annular
  end

  methods (Static)
    function [bsc, weights] = FromMathieu(varargin)
      % Construct a Webber beam using Bessel point-matching
      %
      % Usage
      %   bsc = MathieuBeam(angle, morder, ellipticity, parity, Lmax, ...)
      %   Construct a Webber beam.
      %
      %   [bsc, weights] = WebberBeam(angle, a, parity, Lmax, ...)
      %   Calculate modes and weights for Webber beam.
      %
      % Parameters
      %   - angle (numeric) -- Angle in far-field.
      %
      %   - morder (numeric) -- Mathieu beam mode.
      %   - ellipticity (numeric) -- Ellipticity of Mathieu beam.
      %
      %   - parity (enum) -- Parity of beam (even or odd).
      %
      %   - Lmax (numeric) -- Maximum azimuthal mode number.
      %     If Nmax is specified, defaults to Nmax.
      %     This determines the number of points for point matching.
      %
      % Optional named parameters
      %   - Nmax (numeric) -- Truncation for VSWF coefficients.
      %     For a simpler interface for creating beams without explicit
      %     `Nmax`, see :class:`vswf.PlaneWave`.  Default: ``0``.

      p = inputParser;
      p.addOptional('angle', [], @isnumeric);
      p.addOptional('morder', [], @isnumeric);
      p.addOptional('ellipticity', [], @isnumeric);
      p.addOptional('parity', [], @(x) any(strcmpi(x, {'even', 'odd'})));
      p.addOptional('Lmax', [], @isnumeric);
      p.addParameter('Nmax', 0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      Lmax = p.Results.Lmax;
      if isempty(Lmax)
        Lmax = p.Results.Nmax;
      end

      assert(isnumeric(Lmax) && isscalar(Lmax), ...
          'Lmax must be numeric scalar');

      morder = p.Results.morder;
      assert(isnumeric(morder) && isscalar(morder), ...
          'morder must be numeric scalar');

      ellip = p.Results.ellipticity;
      assert(isnumeric(ellip) && isscalar(ellip), ...
          'ellipticity must be numeric scalar');

      % Generate points for PM
      % Choose points which avoid 0 and pi
      Npts = 2*Lmax+1;
      phi = linspace(0, 2*pi, Npts+2);
      phi = phi(1:end-1) + (phi(2) - phi(1))./2;

      % Calculate MB
      switch p.Results.parity
        case 'odd'
          A = MathieuFunctions(phi, morder, ellip, 'se');
        case 'even'
          A = MathieuFunctions(phi, morder, ellip, 'ce');
        otherwise
          error('Unknown parity value, must be even or odd');
      end

      % Calculate weights
      lmode = -Lmax:Lmax;
      bval = exp(1i.*lmode.*phi.');
      weights = bval \ A.';

      % Calculate Bessel modes
      Etp = ones(2, numel(lmode));
      bsc = ott.bsc.Annular.FromBessel(p.Results.Nmax, ...
          p.Results.angle, Etp, lmode, unmatched{:});

      % Apply weights if needed
      if nargout == 1
        bsc = bsc .* weights.';
      end
    end

    function [bsc, weights] = FromWebber(varargin)
      % Construct a Webber beam using Bessel point-matching
      %
      % Usage
      %   bsc = WebberBeam(angle, a, parity, Lmax, ...)
      %   Construct a Webber beam.
      %
      %   [bsc, weights] = WebberBeam(angle, a, parity, Lmax, ...)
      %   Calculate modes and weights for Webber beam.
      %
      % Parameters
      %   - angle (numeric) -- Angle in far-field.
      %
      %   - a (numeric) -- Parameter describing Webber beam.
      %
      %   - parity (enum) -- Parity of beam (even or odd).
      %
      %   - Lmax (numeric) -- Maximum azimuthal mode number.
      %     If Nmax is specified, defaults to Nmax.
      %     This determines the number of points for point matching.
      %
      % Optional named parameters
      %   - Nmax (numeric) -- Truncation for VSWF coefficients.
      %     For a simpler interface for creating beams without explicit
      %     `Nmax`, see :class:`vswf.PlaneWave`.  Default: ``0``.

      p = inputParser;
      p.addOptional('angle', [], @isnumeric);
      p.addOptional('a', [], @isnumeric);
      p.addOptional('parity', [], @(x) any(strcmpi(x, {'even', 'odd'})));
      p.addOptional('Lmax', [], @isnumeric);
      p.addParameter('Nmax', 0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      Lmax = p.Results.Lmax;
      if isempty(Lmax)
        Lmax = p.Results.Nmax;
      end

      assert(isnumeric(Lmax) && isscalar(Lmax), ...
          'Lmax must be numeric scalar');

      a = p.Results.a;
      assert(isnumeric(a) && isscalar(a), ...
          'a must be numeric scalar');

      % Generate points for PM
      % Choose points which avoid 0 and pi
      Npts = 2*Lmax+1;
      phi = linspace(0, 2*pi, Npts+2);
      phi = phi(1:end-1) + (phi(2) - phi(1))./2;

      % Equation (6) from
      % https://iopscience.iop.org/article/10.1088/1367-2630/14/3/033018
      Ae = 1./(2*sqrt(pi*abs(sin(phi)))) .* exp(1i*a.*log(abs(tan(phi./2))));

      % Add parity
      % TODO: We can probably do this smartly with mmodes
      switch p.Results.parity
        case 'even'
          A = Ae;
        case 'odd'
          A = [Ae(phi <= pi)./1i, -Ae(phi > pi)./1i];
        otherwise
          error('Unknown value for parity');
      end

      % Calculate weights
      lmode = -Lmax:Lmax;
      bval = exp(1i.*lmode.*phi.');
      weights = bval \ A.';

      % Calculate Bessel modes
      Etp = ones(2, numel(lmode));
      bsc = ott.bsc.Annular.FromBessel(p.Results.Nmax, ...
          p.Results.angle, Etp, lmode, unmatched{:});

      % Apply weights if needed
      if nargout == 1
        bsc = bsc .* weights.';
      end
    end

    function [beam, data] = FromBessel(Nmax, theta, Etp, lmode, varargin)
      % Construct a Annular beam from the specified Bessel parameters.
      %
      % Usage
      %   [beam, data] = ott.bsc.Annular.FromBessel(...
      %       Nmax, theta, Etp, lmode, ...)
      %
      % Parameters
      %   - Nmax (numeric) -- Size of beam shape coefficient data.
      %   - theta (N numeric) -- Annular angle in radians from -z direction.
      %   - Etp (2xN numeric) -- Theta and Phi field amplitudes.
      %   - lmode (N numeric) -- Orbital angular momentum number.
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      assert(isnumeric(Nmax) && isscalar(Nmax) ...
          && Nmax >= 0 && round(Nmax) == Nmax, ...
          'Nmax should be an single positive integer');
      assert(isnumeric(Etp) && ismatrix(Etp) && size(Etp, 1) == 2, ...
        'Etp must be 2xN numeric matrix');
        
      % Ensure size of inputs match
      [theta, lmode, Etheta, Ephi] = ott.utils.matchsize(...
          theta(:), lmode(:), Etp(1, :).', Etp(2, :).');

      %% calculate the mode indices we are going to find.
      nTheta = length(theta);
      a = zeros((Nmax*(Nmax+2)), nTheta);
      b = zeros((Nmax*(Nmax+2)), nTheta);

      indexes= (1:nTheta).';

      data = p.Results.data;

      for n = 1:Nmax

        ci_index_m=find(abs(lmode)<=n);
        indt=n+lmode(ci_index_m)+1;

        %power normalisation.
        Nn = 1/sqrt(n*(n+1));

        %Generate the farfield components of the VSWFs
        ci = ott.utils.combined_index(n, -n:n);
        [~, dtY, dpY, data] = data.evaluateYtp(ci, theta, 0);

        %slow indexing.
        szA=sub2ind(size(a),(n-1)*(n+1)+indt,indexes(ci_index_m));
        szY=sub2ind(size(dtY),indt,indexes(ci_index_m));

        %equivalent to dot((1i)^(n+1)*C,E);
        a(szA) = 4*pi*Nn*(-1i)^(n+1) ...
            *(conj(dpY(szY)).*Etheta(ci_index_m) ...
            - conj(dtY(szY)).*Ephi(ci_index_m));
        %equivalent to dot((1i)^(n)*B,E);
        b(szA) = 4*pi*Nn*(-1i)^(n)  ...
            *(conj(dtY(szY)).*Etheta(ci_index_m) ...
            + conj(dpY(szY)).*Ephi(ci_index_m));
      end

      beam = ott.bsc.Annular(a, b, theta);
      beam = beam.makeSparse();
    end
  end

  methods
    function beam = Annular(varargin)
      % Construct a new beam object
      %
      % Usage
      %   beam = Annular() Construct an empty Bsc beam.
      %
      %   beam = Annular(a, b, theta, ...) constructs beam from a/b
      %   coefficients and annular angle vector.
      %   Does not verify that a/b coefficients match theta.
      %   See :meth:`FromBessel` for a Bessel-like constructor.
      %
      % Parameters
      %   a,b (MxN numeric) -- Vectors of VSWF coefficients.
      %   theta (N numeric) -- Annular angles for each beam.

      p = inputParser;
      p.addOptional('a', [], @isnumeric);
      p.addOptional('b', [], @isnumeric);
      p.addOptional('theta', [], @isnumeric);
      p.parse(varargin{:});

      beam = beam@ott.bsc.Bsc(p.Results.a, p.Results.b);
      for ii = 1:numel(beam)
        beam(ii).theta = p.Results.theta(ii);
      end
    end

    function beam = translateZ(beam, z, varargin)
      % Apply translation along z using a phase shift.
      %
      % Usage
      %   beam = beam.translateZ(z)

      ott.utils.nargoutCheck(beam, nargout);
      
      % Make sure sizes match
      Npos = numel(z);
      Nbeams = numel(beam);
      assert(Npos == Nbeams || Npos == 1 || Nbeams == 1, ...
        'Number of beams and positions must match or be scalar');
      if Npos == 1, z = repelem(z, 1, Nbeams); end
      if Nbeams == 1, beam = repmat(beam, 1, Npos); end

      assert(isnumeric(z) && isvector(z), ...
          'z must be numeric vector');
        
      dz = z .* cos([beam.theta]);
      dz = exp(-1i.*dz.*2*pi);

      % Apply translation
      beam = beam .* dz;
    end
  end

  methods % Getters/setters
    function beam = set.theta(beam, val)
      assert(isscalar(val) && isnumeric(val), ...
          'theta must be numeric scalar');
      beam.theta = val;
    end
  end
end
