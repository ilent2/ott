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

  properties
    theta         % Angle describing annular
  end

  methods (Static)
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

      %% calculate the mode indices we are going to find.
      nTheta = length(theta);
      a = zeros((Nmax*(Nmax+2)), nTheta);
      b = zeros((Nmax*(Nmax+2)), nTheta);

      indexes= (1:nTheta).';
      lmode = lmode.';
      theta = theta.';
      Etheta = Etp(1, :).';
      Ephi = Etp(2, :).';

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
      beam.theta = p.Results.theta;
    end

    function beam = translateZ(beam, z, varargin)
      % Apply translation along z using a phase shift.
      %
      % Usage
      %   beam = beam.translateZ(z)

      ott.utils.nargoutCheck(beam, nargout);

      assert(isnumeric(z) && isvector(z), ...
          'z must be numeric vector');

      ibeam = beam;

      for ii = 1:numel(z)
        dz = z(ii) .* cos(ibeam.theta);
        dz = exp(-1i.*dz.*2*pi);

        % Apply translation
        beam(ii) = ibeam .* dz;
      end
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
