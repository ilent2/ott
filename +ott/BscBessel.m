classdef BscBessel < ott.Bsc
%BscBessel representation of a bessel beam and bessel-like beams with OAM
%
% BscBessel properties:
%
% Based on bsc_bessel.m from ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    theta
  end

  methods
    function beam = BscBessel(nmax, theta, varargin)
      %BSCBESSEL construct a new bessel beam or bessel-like beam
      %
      % BSCBESSEL(nmax, theta) creates a new beam with size Nmax
      % for angles theta.
      %
      % BSCBESSEL(..., 'lmode', lmode) specifies the orbital angualr
      % momentum to apply to the beam, default 0.
      %
      % BSCBESSEL(..., 'k_medium', k) specifies the medium wavenumber,
      % defaults to 2*pi, i.e. unit wavelength.

      beam = beam@ott.Bsc();

      % Parse optional inputs
      p = inputParser;
      p.addParameter('polarisation', [ 1 0 ]);
      p.addParameter('lmode', 0);
      p.addParameter('debug', false);
      p.addParameter('k_medium', 2*pi);
      p.parse(varargin{:});

      theta = theta(:);

      beam.type = 'incident';
      beam.basis = 'regular';
      beam.k_medium = p.Results.k_medium;
      beam.theta = theta.';

      Etheta = p.Results.polarisation(:, 1);
      Ephi = p.Results.polarisation(:, 2);
      lmode = p.Results.lmode;
      debug = p.Results.debug;

      szT=size(theta);
      szE=size(Etheta);

      phi=zeros(szT);
      if ~debug
          if nargin==3
              lmode=0;
              if szT(1)~=szE(1)
                  if szE(1)==1
                      Etheta=repmat(Etheta,[szT(1),1]);
                  end
              end
              Ephi=Etheta(:,2);
          end
          
          if nargin==4
              if szE(2)==2
                  if (szE(1)==1)
                      Etheta=repmat(Etheta,[szT(1),1]);
                  end
                  if numel(Ephi)==szT(1)
                      lmode=Ephi;
                  else
                      if numel(Ephi)==1
                          lmode=Ephi(:)*ones(szT(1),1);
                      else
                          error('ott:bsc_bessel:badinput','Input is not valid.')
                      end
                  end
              else
                  lmode=0;
              end
          end
      end

      [nn,mm]=meshgrid([1:nmax],unique(lmode));

      if szE(2)==2
          [nn,mm]=meshgrid([1:nmax],unique([lmode(:)+1;lmode(:)-1]));
          
          m_out_of_bounds=(abs(mm(:))>nn(:));
          nn(m_out_of_bounds)=[];
          mm(m_out_of_bounds)=[];
          
          phi=(phi(:)+1)*[0,pi/2,pi,3*pi/2];
          phi=phi(:);
          
          if numel(lmode)==szT(1)
              lmode=(lmode(:))*[1,1,1,1];
          end
          
          Etheta=repmat(Etheta,[4,1]);
          theta=repmat(theta,[4,1]);
          
          Ephi=(-sin(phi).*Etheta(:,1)+cos(phi) ...
              .*Etheta(:,2)).*exp(1i*lmode(:).*phi(:));
          Etheta=(sign(cos(theta)).*(cos(phi).*Etheta(:,1) ...
              +sin(phi).*Etheta(:,2))).*exp(1i*lmode(:).*phi(:));
      end

      nn=nn(:);
      mm=mm(:);

      Ephi=Ephi(:);
      Etheta=Etheta(:);

      plane = ott.BscPlane(theta, phi, ...
          'polarisation', [Etheta, Ephi], 'Nmax', nmax);
      a = plane.a;
      b = plane.b;

      ci=ott.utils.combined_index(nn,mm);

      a=a(ci,:);
      b=b(ci,:);

      a=reshape(sum(reshape(a,[[szT(1),1/szT(1)].*size(a)]),2) ...
          / size(theta,1)*szT(1),[length(nn),szT(1)]);
      b=reshape(sum(reshape(b,[[szT(1),1/szT(1)].*size(b)]),2) ...
          / size(theta,1)*szT(1),[length(nn),szT(1)]);

      % Make the beam vector and store the coefficients
      [beam.a, beam.b] = ott.Bsc.make_beam_vector(a, b, nn, mm);
    end
  end
end
