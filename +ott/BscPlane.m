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
      %  BSCPLANE(..., 'Nmax', Nmax) and BSCPLANE(..., 'ka', ka)
      %  specify the region where the plane wave beam is valid.
      %
      %  BSCPLANE(..., 'polarisation', [ Etheta Ephi ]) specifies
      %  the polarisation in the theta and phi directions.

      % TODO: The original bsc_plane.m function included support
      % for array inputs theta, phi for multiple beams simulataniously.
      % Is this something we want to support?

      beam = beam@ott.Bsc('incomming', varargin{:});

      % Parse inputs
      p = inputParser;
      p.addParameter('polarisation', [ 1 1i ]);
      p.addParameter('Nmax', []);
      p.addParameter('ka', []);
      p.parse(varargin{:});

      % Store inputs
      beam.theta = theta;
      beam.phi = phi;
      beam.polarisation = p.Results.polarisation;

      % Parse ka/Nmax
      if isempty(p.Results.Nmax) && isempty(p.Results.ka)
        warning('ott:BscPlane:no_range', ...
            'No range specified for plane wave, using range of ka=1');
        Nmax = ott.utils.ka2nmax(1);
      elseif ~isempty(p.Results.ka)
        Nmax = ott.utils.ka2nmax(p.Results.ka);
      elseif ~isempty(p.Results.Nmax)
        Nmax = p.Results.Nmax;
      else
        error('ott:BscPlane:duplicate_range', ...
            'Nmax and ka specified.  Must specify only one');
      end

      ablength = ott.utils.combined_index(Nmax, Nmax);

      a = zeros(ablength, 1);
      b = zeros(ablength, 1);

      for n = 1:Nmax
        iter=[(n-1)*(n+1)+1:n*(n+2)];
        leniter=2*n+1;

        %expand theta and phi components of field to match spherical harmonics
        ET=repmat(Etheta,[1,leniter]);
        EP=repmat(Ephi,[1,leniter]);

        %power normalisation.
        Nn = 1/sqrt(n*(n+1));

        %Generate the farfield components of the VSWFs
        [~,dtY,dpY] = spharm(n,[-n:n],theta,phi);

        %equivalent to dot((1i)^(n+1)*C,E);
        a(iter,:) = 4*pi*Nn*(-1i)^(n+1)*(conj(dpY).*ET - conj(dtY).*EP).';
        %equivalent to dot((1i)^(n)*B,E);
        b(iter,:) = 4*pi*Nn*(-1i)^(n)  *(conj(dtY).*ET + conj(dpY).*EP).';
      end

      p = abs(a).^2 + abs(b).^2;
      non_zeros = p > 1e-15*max(p);

      a(~non_zeros) = 0;
      b(~non_zeros) = 0;

      % Store the coefficients
      beam.a = sparse(a);
      beam.b = sparse(b);
    end

    function beam = translateZ(beam, z)
      %TRANSLATEZ translate a beam along the z-axis

      % Add a warning when the beam is translated outside nmax2ka(Nmax)
      dz = dz + abs(z);
      if dz > ott.utils.nmax2ka(beam.Nmax)/beam.k_medium
        warning('ott:BscPlane:translateZ:outside_nmax', ...
            'Repeated translation of beam outside Nmax region');
      end

      [A, B] = ott.utils.translate_z(beam.Nmax, z);
      beam = blkdiag(A, B) * beam;
    end
  end

end
