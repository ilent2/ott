classdef AxisymShape < ott.shapes.Shape
%AxisymShape abstract class for axisymetric particles
%
% Methods:
%   boundarypoints  calculate boudary points for surface integral
%
% Abstract methods:
%   radii           Calculates the particle radii for angular coordinates
%   normals         Calculates the particle normals for angular coorindates
%   axialSymmetry   Returns x, y, z rotational symmetry (0 for infinite)
%   mirrorSymmetry  Returns x, y, z mirror symmetry
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
  end

  properties (Dependent)
    perimiter
  end

  methods (Abstract)
    radii(shape, theta, phi);
    normals(shape, theta, phi);
    axialSymmetry(shape);
    mirrorSymmetry(shape);

    get_perimiter(shape, varargin);
  end

  methods

    function p = get.perimiter(shape)
      % Get the perimiter of the object
      p = shape.get_perimiter();
    end

    function [rtp, n, ds] = boundarypoints(shape, varargin)
      % BOUNDARYPOINTS calculatse boudnary points for surface integral
      %
      % BOUDNARYPOINTS(npts) or BOUNDARYPOINTS('Nmax', Nmax) calculates
      % the boundary points specifying the number of points to calculate
      % or the Nmax to calculate points for.

      % Parse inputs
      p = inputParser;
      p.addOptional('npts', []);
      p.addParameter('Nmax', []);
      p.parse(varargin{:});

      % Determine the number of points to use
      if isempty(p.Results.npts) && isempty(p.Results.Nmax)
        error('Must specify either npts or Nmax');
      elseif ~isempty(p.Results.npts) && ~isempty(p.Results.Nmax)
        error('Both number of points and Nmax specified');
      elseif ~isempty(p.Results.npts)
        ntheta = p.Results.npts;
      else
        Nmax = p.Results.Nmax;
        ntheta = ceil(ott.utils.combined_index(Nmax, Nmax).^2/20+5);
      end

      % Check that the point is axissymetric
      axisym = shape.axialSymmetry();
      if axisym(3) ~= 0
        error('Only supports axisymetric particles');
      end

      % Calculate the point spacing
      ds = shape.perimiter / ntheta / 2.0;
      s = 0:ds:(shape.perimiter/2.0);

      % TODO: This isn't optimal for all shapes
      % We should be able to optimise this interface
      % TODO: This won't work if we don't provide Nmax
      % Why doesn't this work with more points?
      [theta, phi] = ott.utils.angulargrid(4*(Nmax + 2), 1);
      xyz = shape.locations(theta, phi);
      rho = xyz(:, 1);
      z = xyz(:, 3);

      % The following is based on axisym_boundarypoints from OTTv1

      %Don't ask me how this works. It does. It's simple algebra in the end...
      %and yes, it can be done more intelligently.
      zout=zeros(ntheta,1);
      rhoout=zout;
      nxyz=zeros(ntheta,3);

      sdeficit=0;
      ncum=0;

      for ii=2:length(rho)
          N=s(ii-1)/ds;
          Nused=round(N+sdeficit);

          nc=[-(z(ii)-z(ii-1)),0,(rho(ii)-rho(ii-1))];
          nc=nc/norm(nc,2);

          if Nused>=1
              drho=(rho(ii)-rho(ii-1))/N*ones(Nused,1);

              rhot=cumsum(drho)-drho/2-sdeficit*drho(1);
              rhoout([ncum+[1:Nused]])=rho(ii-1)+rhot;

              dz=(z(ii)-z(ii-1))/N*ones(Nused,1);

              zt=cumsum(dz)-dz/2-sdeficit*dz(1);
              zout([ncum+[1:Nused]])=z(ii-1)+zt;

              nxyz([ncum+[1:Nused]],:)=repmat(nc,[length(zt),1]);

              sdeficit=(N-Nused+sdeficit);
          else
              sdeficit=sdeficit+N;
          end

          ncum=ncum+Nused;

      end

      %converts the cylindrical coordinates into spherical coordinates
      [n,rtp]=ott.utils.xyzv2rtpv(nxyz,[rhoout,zeros(size(rhoout)),zout]);

      dst=zeros(size(rtp, 1),3);

      %calcultes area elements
      dst(2:end-1,1)=(rhoout(3:end)-rhoout(1:end-2))/2;
      dst(2:end-1,3)=(zout(3:end)-zout(1:end-2))/2;
      dst(1,1)=(mean(rhoout(1:2))-rho(1));
      dst(1,3)=(mean(zout(1:2))-z(1));
      dst(end,1)=(rho(end)-mean(rhoout(end-1:end)));
      dst(end,3)=(z(end)-mean(zout(end-1:end)));

      % a general axisymmetric conic region has the
      % following area (sans factor of 2*pi):
      ds=(rtp(:,1).*sqrt(sum(abs(dst).^2,2)).*sin(rtp(:,2)));
    end
  end
end
