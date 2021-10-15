classdef BscPlane < ott.Bsc
%BscPlane representation of a plane wave in VSWF coefficients
%
% BscPlane properties:
%   theta           Beam direction (polar angle)
%   phi             Beam direction (azimuthal angle)
%   polarisation    Beam polarisation [ Etheta Ephi ]
%
% BscPlane methods:
%   translateZ      Translates the beam and checks within beam range
%
% Based on bsc_plane.m from ottv1.
%
% See also BscPlane and ott.Bsc.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    theta           % Beam direction (polar angle)
    phi             % Beam direction (azimuthal angle)
    polarisation    % Beam polarisation [ Etheta Ephi ]
  end

  methods
    function beam = BscPlane(theta, phi, varargin)
      %BSCPLANE construct a new plane wave beam
      %
      %  BSCPLANE(theta, phi) creates a new circularly polarised plane
      %  wave beam in a direction specified by theta and phi.
      %     theta polar angle from +z axis (rad)
      %     phi   azimuthal angle, measured from +x towards +y axes (rad)
      %
      % theta and phi can be arrays, in which case multiple VSWF
      % expansions are calculated for each angle.
      %
      %  BSCPLANE(..., 'k_medium', k) specifies the wavenumber in the
      %  medium.  Defaults to 2*pi, i.e. wavelength = 1.
      %
      %  BSCPLANE(..., 'Nmax', Nmax) and BSCPLANE(..., 'radius', a)
      %  specify the region where the plane wave beam is valid.
      %
      %  BSCPLANE(..., 'polarisation', [ Etheta Ephi ]) specifies
      %  the polarisation in the theta and phi directions.

      beam = beam@ott.Bsc();

      % Parse inputs
      p = inputParser;
      p.addParameter('polarisation', [ 1 1i ]);
      p.addParameter('Nmax', []);
      p.addParameter('radius', []);
      p.addParameter('wavelength0', 1);
      p.addParameter('power', []);

      p.addParameter('omega', 2*pi);
      p.addParameter('k_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength_medium', []);

      p.parse(varargin{:});

      beam.type = 'incident';
      beam.basis = 'regular';
      beam.k_medium = ott.Bsc.parser_k_medium(p, 2*pi);
      beam.omega = p.Results.omega;

      % If points aren't specified explicitly, use meshgrid
      if length(theta) ~= length(phi)
        [theta, phi] = meshgrid(theta, phi);
        theta = theta(:);
        phi = phi(:);
      end

      % Store inputs
      beam.theta = theta;
      beam.phi = phi;

      % Store polarisation
      beam.polarisation = p.Results.polarisation;
      if size(p.Results.polarisation, 1) ~= length(phi) ...
          && size(p.Results.polarisation, 1) ~= 1
        error('Polarisation must be either 1x2 or Nx2 (N = # beams)');
      end
      if size(p.Results.polarisation, 1) ~= length(phi)
        beam.polarisation = repmat(beam.polarisation, length(phi), 1);
      end

      % Parse ka/Nmax
      if isempty(p.Results.Nmax) && isempty(p.Results.radius)
        warning('ott:BscPlane:no_range', ...
            'No range specified for plane wave, using range of ka=1');
        Nmax = ott.utils.ka2nmax(1);
      elseif ~isempty(p.Results.radius)
        Nmax = ott.utils.ka2nmax(p.Results.radius * beam.k_medium);
      elseif ~isempty(p.Results.Nmax)
        Nmax = p.Results.Nmax;
        assert(isscalar(Nmax), 'Nmax must be scalar');
      else
        error('ott:BscPlane:duplicate_range', ...
            'Nmax and ka specified.  Must specify only one');
      end

      ablength = ott.utils.combined_index(Nmax, Nmax);

      a = zeros(ablength, length(beam.theta));
      b = zeros(ablength, length(beam.theta));

      Etheta = beam.polarisation(:, 1);
      Ephi = beam.polarisation(:, 2);

      for n = 1:Nmax
        iter=[(n-1)*(n+1)+1:n*(n+2)];
        leniter=2*n+1;

        %expand theta and phi components of field to match spherical harmonics
        ET=repmat(Etheta,[1,leniter]);
        EP=repmat(Ephi,[1,leniter]);

        %power normalisation.
        Nn = 1/sqrt(n*(n+1));

        %Generate the farfield components of the VSWFs
        [~,dtY,dpY] = ott.utils.spharm(n,[-n:n],theta,phi);

        %equivalent to dot((1i)^(n+1)*C,E);
        a(iter,:) = 4*pi*Nn*(-1i)^(n+1)*(conj(dpY).*ET - conj(dtY).*EP).';
        %equivalent to dot((1i)^(n)*B,E);
        b(iter,:) = 4*pi*Nn*(-1i)^(n)  *(conj(dtY).*ET + conj(dpY).*EP).';
      end

      pw = abs(a).^2 + abs(b).^2;
      non_zeros = pw > 1e-15*max(pw);

      a(~non_zeros) = 0;
      b(~non_zeros) = 0;

      % Store the coefficients
      beam.a = sparse(a);
      beam.b = sparse(b);

      % Normalize the beam power
      if ~isempty(p.Results.power)
        beam.power = p.Results.power;
      end
    end

    function [beam, A, B] = translateZ(beam, varargin)
      %TRANSLATEZ translate a beam along the z-axis
      %
      % TRANSLATEZ(z) translates by a distance z along z axis.
      %
      % [beam, A, B] = TRANSLATEZ(z) returns the translation matrices
      % and translated beam.
      %
      % [beam, AB] = TRANSLATEZ(z) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam
      
      % TODO: This could be implemented as a phase shift, or should
      % it be implementedin translateXyz instead?

      p = inputParser;
      p.addOptional('z', []);
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      if ~isempty(p.Results.z)
        z = p.Results.z;

        % Add a warning when the beam is translated outside nmax2ka(Nmax)
        beam.dz = beam.dz + abs(z);
        if beam.dz > ott.utils.nmax2ka(beam.Nmax)/beam.k_medium
          warning('ott:BscPlane:translateZ:outside_nmax', ...
              'Repeated translation of beam outside Nmax region');
        end

        % Convert to beam units
        z = z * beam.k_medium / 2 / pi;

        [A, B] = beam.translateZ_type_helper(z, [p.Results.Nmax, beam.Nmax]);

      else
        error('Wrong number of arguments');
      end

      % Apply the translation
      beam = beam.translate(A, B);
      beam.basis = 'regular';

      % Pack the rotated matricies into a single ABBA object
      if nargout == 2
        A = [ A B ; B A ];
      end
    end
    
    function [sbeam, beam] = scatter(beam, tmatrix, varargin)
      %SCATTER scatter a beam using a T-matrix
      %
      % For BscPlane: Adds an extra check to make sure the beam Nmax is
      % large enough to represent a plane wave.
      %
      % For full documentation see Bsc/scatter.

      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax2 = 0;
      for ii = 1:numel(tmatrix)
        maxNmax2 = max(maxNmax2, tmatrix(ii).Nmax(2));
      end
      
      if beam.Nmax < maxNmax2
        warning('ott:BscPlane:scatter:small_nmax', ...
          'The beam Nmax may be too small for the current T-matrix');
      end
      
      % Do the default processing
      [sbeam, beam] = scatter@ott.Bsc(beam, tmatrix, varargin{:});
    end
  end
end
