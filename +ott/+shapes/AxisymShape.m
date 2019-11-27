classdef AxisymShape < ott.shapes.Shape
%AxisymShape abstract class for axisymetric particles
%
% Methods
%   - boundarypoints  calculate boudary points for surface integral
%
% Abstract methods
%   - radii           Calculates the particle radii for angular coordinates
%   - normals         Calculates the particle normals for angular coorindates
%   - axialSymmetry   Returns x, y, z rotational symmetry (0 for infinite)

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

    get_perimiter(shape, varargin);
  end

  methods (Access=protected)

    function npts = boundarypoints_npts(shape, varargin)
      % Helper for boundarypoints

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
        npts = p.Results.npts;
      else
        Nmax = p.Results.Nmax;
        npts = ceil(ott.utils.combined_index(Nmax, Nmax).^2/20+5);
      end

    end

    function ds = boundarypoints_area(shape, rho, z, rhoout, zout, rtp)
      % Helper for boundarypoints

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

    function [rtp, n, ds] = boundarypoints_rhoz(shape, rho, z, varargin)
      % Helper for boundarypoints

      ntheta = shape.boundarypoints_npts(varargin{:});

      % Check that the point is axissymetric
      axisym = shape.axialSymmetry();
      if axisym(3) ~= 0
        error('Only supports axisymetric particles');
      end

      % Calculate the point spacing
      ds = shape.perimiter / ntheta / 2.0;

      % The following is based on axisym_boundarypoints from OTTv1

      % Calculate length of each line segment
      s=sqrt((rho(2:end)-rho(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2);

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
              rhoout(ncum+(1:Nused))=rho(ii-1)+rhot;

              dz=(z(ii)-z(ii-1))/N*ones(Nused,1);

              zt=cumsum(dz)-dz/2-sdeficit*dz(1);
              zout(ncum+(1:Nused))=z(ii-1)+zt;

              nxyz(ncum+(1:Nused),:)=repmat(nc,[length(zt),1]);

              sdeficit=(N-Nused+sdeficit);
          else
              sdeficit=sdeficit+N;
          end

          ncum=ncum+Nused;

      end
      
      % Truncate points if not allocated
      if ncum < ntheta
        warning('OTT:SHAPES:AXISYMSHAPE:boundarypoints_length', ...
          'Number of points generated does not match request');
        zout = zout(1:ncum);
        rhoout = zout(1:ncum);
        nxyz = nxyz(1:ncum, :);
      end

      %converts the cylindrical coordinates into spherical coordinates
      [n,rtp]=ott.utils.xyzv2rtpv(nxyz,[rhoout,zeros(size(rhoout)),zout]);

      % Calculate area elements
      ds = shape.boundarypoints_area(rho, z, rhoout, zout, rtp);

    end
  end

  methods

    function p = get.perimiter(shape)
      % Get the perimiter of the object
      p = shape.get_perimiter();
    end

    function [rtp, n, ds] = boundarypoints(shape, varargin)
      % BOUNDARYPOINTS calculates boundary points for surface integral
      %
      % [rtp, n, ds] = BOUDNARYPOINTS(npts) calculates the boundary points
      % and surface normal vectors in spherical coordinates and the area
      % elements of each ring.
      %
      % BOUNDARYPOINTS('Nmax', Nmax) takes a guess at a suitable npts
      % for the given Nmax.

      npts = shape.boundarypoints_npts(varargin{:});

      % Sample points along the boundary, this isn't needed if
      % we have a shape described by a finite set of boundary points
      % TODO: How many points should we sample?
      % Angular grid doesn't include top and bottom, so we add them
      [theta, phi] = ott.utils.angulargrid(npts*2, 1);
      theta = [0.0; theta; pi];
      phi = [phi(1); phi; phi(end)];

      xyz = shape.locations(theta, phi);
      rho = xyz(:, 1);
      z = xyz(:, 3);

      [rtp, n, ds] = shape.boundarypoints_rhoz(rho, z, varargin{:});
    end
  end
end
