classdef BscBessel < ott.Bsc
%BscBessel representation of a Bessel beam and Bessel-like beams with OAM
%
% BscBessel properties:
%   theta        Bessel beam angle in far-field
%   a            (Bsc) Beam shape coefficients a vector
%   b            (Bsc) Beam shape coefficients b vector
%   type         (Bsc) Beam type (incident, scattered, total)
%   basis        (Bsc) VSWF beam basis (incoming, outgoing or regular)
%   Nmax         (Bsc) Truncation number for VSWF coefficients
%   power        (Bsc) Power of the beam [M*L^2/S^3]
%   Nbeams       (Bsc) Number of beams in this Bsc object
%   wavelength   (Bsc) Wavelength of beam [L]
%   speed        (Bsc) Speed of beam in medium [L/T]
%   omega        (Bsc) Angular frequency of beam [2*pi/T]
%   k_medium     (Bsc) Wavenumber in medium [2*pi/L]
%   dz           (Bsc) Absolute cumulative distance the beam has moved
%
% See also BscBessel, ott.Bsc and ott.BscPlane.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    theta       % Bessel beam angle in far-field
    lmode       % Orbital angular momentum mode number
  end

  methods
    function beam = BscBessel(nmax, theta, varargin)
      %BSCBESSEL construct a new Bessel beam or Bessel-like beam
      %
      % BSCBESSEL(Nmax, theta, ...) creates a new beam with size Nmax
      % for angles theta (in radians).
      %
      % If theta is an array, creates a Bsc object with multiple
      % beams, see Nbeams and mergeBeams.  Most of the optional
      % parameters support vector or matrix input with Nbeams rows.
      %
      % Optional named arguments:
      %   lmode    num    orbital angular momentum, default 0
      %   k_medium num    Wavenumber in medium, default 2*pi
      %   polarisation [x, y]  beam polarisation, deafult [1, 0]
      %   Etheta   num    beam polarisation (spherical coords)
      %   Ephi     num    beam polarisation (spherical coords)
      %
      % Only polarisation or both Etheta and Ephi should be specified.
      %
      % See also beam, mergeBeams, ott.BscPlane

      beam = beam@ott.Bsc();

      % Parse optional inputs
      p = inputParser;
      p.addParameter('polarisation', []);
      p.addParameter('Etheta', []);
      p.addParameter('Ephi', []);
      p.addParameter('lmode', 0);
      p.addParameter('k_medium', 2*pi);
      p.parse(varargin{:});

      % Ensure theta is a column vector
      theta = theta(:);

      % Check theta
      if any(abs(theta) < 1e-6)
        warning('OTT:BscBessel:theta_zero', ...
          'theta = 0 is not a bessel beam, may have unexpected results');
      end

      % Reshape lmode (dupliate bellow if needed)
      lmode = p.Results.lmode(:);

      % Get Etheta and Ephi and re-map to bsc_bessel argument names
      if isempty(p.Results.polarisation) ...
          && ~isempty(p.Results.Etheta) && ~isempty(p.Results.Ephi)
        Etheta = p.Results.Etheta(:);
        Ephi = p.Results.Ephi(:);

        try
          [theta, lmode, Etheta, Ephi] = ott.utils.matchsize(...
              theta, lmode, Etheta, Ephi);
        catch ME
          error(['Number of elements in theta, lmode, Etheta, ' ...
              'Ephi should be equal or 1'], 'OTT:BscBessel:size_mismatch');
        end

        indexes=[1:length(theta)].';

      elseif isempty(p.Results.Etheta) && isempty(p.Results.Ephi)

        % Default polarisation
        polarisation = p.Results.polarisation;
        if isempty(polarisation)
          polarisation = [1, 0];
        end

        assert(size(polarisation, 2) == 2, ...
            'polarisation must be 2 column matrix or vector');

        try
          [theta, lmode, polarisation] = ott.utils.matchsize(...
              theta, lmode, polarisation);
        catch ME
          error(['Number of elements in theta, lmode, polarisation ' ...
              'should be equal or 1'], 'OTT:BscBessel:size_mismatch');
        end

        lmode=lmode+[-1,1]; %left is -1i for Ephi, right is +1i for Ephi;
        lmode=lmode(:);
        indexes=[[1:size(polarisation,1)].'; ...
                 [1:size(polarisation,1)].']; %packing for theta uniform

        %need to record polarisation and modify Etheta and Ephi.
        polarisation_weights=([1,1i;1,-1i]*(polarisation.')).'/2;
        Ephi=[-1i*ones(size(polarisation,1),1); ...
               1i*ones(size(polarisation,1),1)].*polarisation_weights(:);
        Etheta=[ones(size(polarisation,1),1); ...
                ones(size(polarisation,1),1)] ...
                .*polarisation_weights(:).*sign(cos(theta(indexes)));
      else
        error('OTT:BscBessel:too_many_polarisations', ...
            'Only polarisation or both Etheta and Ephi should be supplied');
      end

      % preamble
      nTheta=length(theta);

      %% calculate the mode indices we are going to find.
      % [nn,mm]=combined_index([1:nmax*(nmax+2)].');
      a = zeros((nmax*(nmax+2)),nTheta);
      b = zeros((nmax*(nmax+2)),nTheta);

      for n = 1:nmax

        ci_index_m=find(abs(lmode)<=n);
        indt=n+lmode(ci_index_m)+1;

        %power normalisation.
        Nn = 1/sqrt(n*(n+1));

        %Generate the farfield components of the VSWFs
        [~,dtY,dpY] = ott.utils.spharm(n,theta, 0);

        %slow indexing.
        szA=sub2ind(size(a),(n-1)*(n+1)+indt,indexes(ci_index_m));
        szY=sub2ind(size(dtY),indexes(ci_index_m),indt);

        dtY=dtY(:);
        dpY=dpY(:);

        %equivalent to dot((1i)^(n+1)*C,E);
        a(szA) = 4*pi*Nn*(-1i)^(n+1) ...
            *(conj(dpY(szY)).*Etheta(ci_index_m) ...
            - conj(dtY(szY)).*Ephi(ci_index_m));
        %equivalent to dot((1i)^(n)*B,E);
        b(szA) = 4*pi*Nn*(-1i)^(n)  ...
            *(conj(dtY(szY)).*Etheta(ci_index_m) ...
            + conj(dpY(szY)).*Ephi(ci_index_m));

      end

      % Setup the beam object
      beam.type = 'incident';
      beam.basis = 'regular';
      beam.k_medium = p.Results.k_medium;
      beam.theta = theta.';
      beam.lmode = p.Results.lmode.';

      ci=find(any(abs(a)|abs(b),2));
      [nn,mm]=ott.utils.combined_index(ci);
      a=a(ci,:);
      b=b(ci,:);

      % Make the beam vector and store the coefficients
      [beam.a, beam.b] = ott.Bsc.make_beam_vector(a, b, nn, mm);
    end
    
    function [sbeam, beam] = scatter(beam, tmatrix, varargin)
      %SCATTER scatter a beam using a T-matrix
      %
      % For BscBessel: Adds an extra check to make sure the beam Nmax is
      % large enough to represent a bessel beam.
      %
      % For full documentation see Bsc/scatter.

      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax2 = 0;
      for ii = 1:numel(tmatrix)
        maxNmax2 = max(maxNmax2, tmatrix(ii).Nmax(2));
      end
      
      if beam.Nmax < maxNmax2
        warning('ott:BscBessel:scatter:small_nmax', ...
          'The beam Nmax may be too small for the current T-matrix');
      end
      
      % Do the default processing
      [sbeam, beam] = scatter@ott.Bsc(beam, tmatrix, varargin{:});
    end
  end
end
