classdef Ebcm < ott.Tmatrix ...
    & ott.scat.utils.RelativeMediumProperty
% Constructs a T-matrix using extended boundary conditions method.
% Inherits from :class:`ott.Tmatrix`.
%
% Implements the extended boundary conditions methods for rotationally
% symmetric homogeneous particles.
%
% Properties
%   - points         -- (2xN numeric) Surface points [r; theta]
%   - normals        -- (2xN numeric) Normals at points [nr; ntheta]
%   - areas          -- (1xN numeric) Conic section surface areas
%
% Static methosd
%   - FromStarShape   -- Construct from a star shaped object
%   - FromAxisymInterpShape -- Construct using AxisymInterp shape
%
% Additional methods/properties inherited from :class:`Tmatrix`.
%
% This class is based on tmatrix_ebcm_axisym.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    points         % (2xN numeric) Surface points [r; theta]
    normals        % (2xN numeric) Normals at points [nr; ntheta]
    areas          % (1xN numeric) Conic section surface areas
  end

  methods (Static)
    function varargout = FromAxisymInterpShape(varargin)
      % Construct a T-amtrix using EBCM from a AxisymInterp shape object.
      %
      % Usage
      %   tmatrix = Ebcm.FromAxisymInterpShape(shape, relativeMedium, ...)
      %   Calculate external T-matrix (unless `internal` is true)
      %
      %   [external, internal] = Ebcm.FromAxisymInterpShape(...)
      %   Calculate both the internal and external T-matrices.
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- Description of shape geometry.
      %     Object must be a AxisymInterp or be castable to AxisymInterp.
      %
      %   - relativeMedium (ott.beam.medium.RelativeMedium) -- The relative
      %     medium describing the particle's material.
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Used to convert `shape` input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., `shape` is already in relative units).
      %
      % Additional parameters are passed to :meth:`Ebcm`.

      p = inputParser;
      p.addRequired('shape', @(x) isa(x, 'ott.shapes.Shape'));
      p.addRequired('relativeMedium', ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addParameter('wavelength', 1.0);
      p.addParameter('Nmax', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Ensure shape is an AxisymInterp
      shape = p.Results.shape;
      if ~isa(shape, 'ott.shapes.AxisymInterp')
        shape = ott.shapes.AxisymInterp(shape);
      end

      % Rescale shape
      shape = shape ./ p.Results.wavelength;

      % Get Nmax for npts calculation
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        Nmax = ott.utils.ka2nmax(2*pi*shape.maxRadius);
      else
        assert(isnumeric(Nmax) && isscalar(Nmax) && Nmax > 0, ...
            'Nmax must be positive numeric scalar');
      end

      % Calculate desired number of boundary points
      npts = ceil(ott.utils.combined_index(Nmax, Nmax).^2/20+5);

      % Calculate boundary points and change coordinates
      [xyz, nxyz, ds] = shape.boundaryPoints(npts);
      [nrtp, rtp] = ott.utils.xyzv2rtpv(nxyz, xyz);

      % Construct T-matrices
      [varargout{1:nargout}] = ott.scat.vswf.Ebcm(rtp, nrtp, ds, ...
          p.Results.relativeMedium, 'Nmax', Nmax, unmatched{:});
    end

    function varargout = FromStarShape(varargin)
      % Construct a T-matrix using EBCM from a star-shaped object
      %
      % Usage
      %   tmatrix = Ebcm.FromStarShape(shape, relativeMedium, ...)
      %   Calculate external T-matrix (unless `internal` is true)
      %
      %   [external, internal] = Ebcm.FromStarShape(shape, relativeMedium, ...)
      %   Calculate both the internal and external T-matrices.
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- Description of shape geometry.
      %     Object should implement a valid starRadii method.
      %
      %   - relativeMedium (ott.beam.medium.RelativeMedium) -- The relative
      %     medium describing the particle's material.
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Used to convert `shape` input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., `shape` is already in relative units).
      %
      %   - angles (N numeric) -- Angles describing edge of line segments.
      %     The function calls Ebcm with `N-1` points positioned at the
      %     centre of these line segments.
      %
      % All other parameters are passed to the class constructor.

      p = inputParser;
      p.addRequired('shape', @isnumeric);
      p.addRequired('relativeMedium', ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addParameter('wavelength', 1.0);
      p.addParameter('Nmax', []);
      p.addParameter('angles', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Rescale shape
      shape = shape ./ p.Results.wavelength;

      assert(shape.zRotSymmetry == 0, ...
          'shape must be rotationally symmetric about the z axis');

      % Get Nmax for npts calculation
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        Nmax = ott.utils.ka2nmax(2*pi*shape.maxRadius);
      else
        assert(isnumeric(Nmax) && isscalar(Nmax) && Nmax > 0, ...
            'Nmax must be positive numeric scalar');
      end

      angles = p.Results.angles;
      if isempty(angles)
        % Calculate desired number of boundary points
        npts = ceil(ott.utils.combined_index(Nmax, Nmax).^2/20+5);

        % Calculate angles
        angles = linspace(0, pi, npts+1);
      end

      % Calculate radii at points
      phi = zeros(size(angles)));
      radii = shape.starRadii(angles, phi);
      rtp = [radii(:), angles(:), phi(:)].';
      xyz = ott.utils.rtp2xyz(rtp);

      % Calculate normals
      nxyz = shape.normals(xyz);
      nrtp = ott.utils.xyzv2rtpv(nxyz, xyz);

      % Calculate area of segments
      ds = 

      % Construct T-matrices
      [varargout{1:nargout}] = ott.scat.vswf.Ebcm(rtp, nrtp, ds, ...
          p.Results.relativeMedium, unmatched{:});
    end
  end

  methods
    function [tmatrix, internal] = Ebcm(varargin)
      % Calculates T-matrix using extended boundary condition method
      %
      % Usage
      %   tmatrix = Ebcm(points, normals, area, relativeMedium, ...)
      %   Calculate external T-matrix (unless `internal` is true)
      %
      %   [external, internal] = Ebcm(...)
      %   Calculate both the internal and external T-matrices.
      %
      % Parameters
      %   - points (2xN numeric) -- Coordinates describing surface.
      %     Spherical coordinates (omitting azimuthal angle: `[r; theta]`).
      %
      %   - normals (2xN numeric) -- Normals at surface points.
      %
      %   - areas (N numeric) -- Area elements at surface points.
      %
      %   - relativeMedium (ott.beam.medium.RelativeMedium) -- The relative
      %     medium describing the particle's material.
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Used to convert radius input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., radius is already in relative units).
      %
      %   - internal (logical) -- If true, the returned T-matrix is
      %     an internal T-matrix.  Ignored for two outputs.
      %     Default: ``false``.
      %
      %   - Nmax (numeric) -- Size of the VSWF expansion used for the
      %     T-matrix calculation.  In some cases it can be reduced
      %     after construction.
      %     Default: ``ott.utis.ka2nmax(2*pi*shape.maxRadius)`` (may
      %     need different values to give convergence for some shapes).
      %
      %   - verbose (logical) -- If true, outputs the condition number
      %     of the Q and RgQ matrices.

      p = inputParser;
      p.addRequired('points', @isnumeric);
      p.addRequired('normals', @isnumeric);
      p.addRequired('areas', @isnumeric);
      p.addParameter('wavelength', 1.0);
      p.addParameter('internal', false);
      p.addParameter('Nmax', []);
      p.addParameter('verbose', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct base class and store parameters (changing units)
      tmatrix = tmatrix@ott.scat.vswf.Tmatrix(unmatched{:});
      tmatrix.points = p.Results.points;
      tmatrix.points(1, :) = tmatrix.points(1, :) ./ p.Results.wavelength;
      tmatrix.normals = p.Results.normals;
      tmatrix.areas = p.Results.areas ./ p.Results.wavelength.^2;
      tmatrix.relativeMedium = p.Results.relativeMedium;

      % TODO: Continue from here down






      % Parse inputs
      pa = inputParser;
      pa.addParameter('z_mirror_symmetry', false);
      pa.addParameter('internal', false);
      pa.addParameter('invmethod', []);
      pa.addParameter('verbose', false);
      pa.parse(varargin{:});

      % Store inputs k_medium and k_particle
      [tmatrix.k_medium, tmatrix.k_particle] = ...
          tmatrix.parser_wavenumber(pa, 2*pi);

      % Get or estimate Nmax from the inputs
      if isempty(pa.Results.Nmax)
        maxRadius = max(rtp(:, 1));
        if p.Results.internal
          Nmax = ott.utils.ka2nmax(maxRadius * abs(tmatrix.k_particle));
        else
          Nmax = ott.utils.ka2nmax(maxRadius * abs(tmatrix.k_medium));
        end
      else
        Nmax = pa.Results.Nmax;
      end

      % Choose default value for inverse method
      invmethod = pa.Results.invmethod;
      if isempty(invmethod)
        if imag(tmatrix.k_particle) == 0.0
          % TODO: We should also check if k_particle is small
          invmethod = 'inv';
        else
          invmethod = 'forwardslash';
        end
      end

      mirrorsym = pa.Results.z_mirror_symmetry;

      %%%% set up containters
      Y=cell(Nmax,1);
      Ytheta=Y;
      Yphi=Y;

      kr=tmatrix.k_particle*rtp(:,1);
      ikr=1./kr;
      kr_=tmatrix.k_medium*rtp(:,1);
      ikr_=1./kr_;

      jkr=cell(Nmax,1);
      jkr_=jkr;
      hkr_=jkr_;

      djkr=jkr;
      djkr_=djkr;
      dhkr_=djkr_;

      Nn=1./sqrt((1:Nmax).*((1:Nmax)+1));

      Nm=repmat(Nn,[Nmax,1]);
      Nm=Nm.'.*Nm;
      
      % Check kr
      if any(kr == 0)
        error('kr must not be zero');
      end
      
      import ott.utils.*;

      %make all the spharms and bessel functions NOW!!!!
      for ii=1:Nmax

        [Y{ii},Ytheta{ii},Yphi{ii}] = spharm(ii,rtp(:,2),rtp(:,3));

        jkr{ii}=sbesselj(ii,kr);
        jkr_{ii}=sbesselj(ii,kr_);
        hkr_{ii}=sbesselh1(ii,kr_); %%this is conj

        djkr{ii}=sbesselj(ii-1,kr)-ii*jkr{ii}./kr;
        djkr_{ii}=sbesselj(ii-1,kr_)-ii*jkr_{ii}./kr_;
        dhkr_{ii}=sbesselh1(ii-1,kr_)-ii*hkr_{ii}./kr_; %%this is conj

      end

      lengthJs=Nmax^2*(Nmax+2)-sum((1:Nmax-1).*((1:Nmax-1)+1));
      J11=zeros(lengthJs,1);
      J12=J11;
      J21=J11;
      J22=J11;

      if ~pa.Results.internal
        RgJ11=J11;
        RgJ12=J11;
        RgJ21=J11;
        RgJ22=J11;
      end

      i_indexes=J11;
      j_indexes=J11;

      Qv=zeros(lengthJs*4,1);
      RgQv=zeros(lengthJs*4,1);

      jcount=0;

      fillnum=(2*Nmax+1);
      rM=repmat(normals(:,1),[1,fillnum]);
      thetaM=repmat(normals(:,2),[1,fillnum]);
      iKr=repmat(ikr,[1,fillnum]);
      iKr_=repmat(ikr_,[1,fillnum]);
      dS=repmat(ds,[1,fillnum]);

      pi2=2*pi;
      %%%% end setup containers
      
      k_medium = tmatrix.k_medium;
      k_particle = tmatrix.k_particle;

      for jj=1:Nmax
        %remember that when phi = 0, conj(Yphi)=-Yphi, conj(Ytheta)=Ytheta.

        for kk=1:Nmax
          p=min(kk,jj);

          indexzk=kk+1+(-p:p);
          indexzj=jj+1+(-p:p);

          els= 1:length(indexzk);

          jindx=jcount+els;

          i_indexes(jindx)=kk*(kk+1)+(-p:p);
          j_indexes(jindx)=jj*(jj+1)+(-p:p);

          %Calculate the cross products of the spherical harmonics
          YtYp=-Ytheta{jj}(:,indexzj).*Yphi{kk}(:,indexzk);
          YpYt=Yphi{jj}(:,indexzj).*Ytheta{kk}(:,indexzk);

          YYt=Y{jj}(:,indexzj).*Ytheta{kk}(:,indexzk);
          YtY=Ytheta{jj}(:,indexzj).*Y{kk}(:,indexzk);

          YYp=-Y{jj}(:,indexzj).*Yphi{kk}(:,indexzk);
          YpY=Yphi{jj}(:,indexzj).*Y{kk}(:,indexzk);

          YpYp=-Yphi{jj}(:,indexzj).*Yphi{kk}(:,indexzk);
          YtYt=Ytheta{jj}(:,indexzj).*Ytheta{kk}(:,indexzk);

          %calculate the cross products of spherical bessel functions.
          jh_=repmat(jkr{jj}.*(hkr_{kk}),[1,fillnum]);
          if ~pa.Results.internal
            jj_=repmat(jkr{jj}.*jkr_{kk},[1,fillnum]);
          end

          jdh_=repmat(jkr{jj}.*(dhkr_{kk}),[1,fillnum]);
          djh_=repmat(djkr{jj}.*(hkr_{kk}),[1,fillnum]);

          if ~pa.Results.internal
            jdj_=repmat(jkr{jj}.*djkr_{kk},[1,fillnum]);
            djj_=repmat(djkr{jj}.*jkr_{kk},[1,fillnum]);
          end

          djdh_=repmat(djkr{jj}.*(dhkr_{kk}),[1,fillnum]);
          if ~pa.Results.internal
            djdj_=repmat(djkr{jj}.*djkr_{kk},[1,fillnum]);
          end

          % perform the cross product followed by dotting the
          % normal vector and summing over the area elements.

          mfac=1;
          if ~(mirrorsym&&(rem(jj,2)==rem(kk,2)))
            J11(jindx) = mfac.*(pi2*Nm(jj,kk)* ...
                sum(dS(:,els).*rM(:,els).*jh_(:,els).*(YtYp-YpYt),1));
            J22(jindx) = mfac.*(pi2*Nm(jj,kk)* ...
                sum(dS(:,els).*(rM(:,els).*djdh_(:,els).*(YtYp-YpYt) ...
                -thetaM(:,els).*(jj*(jj+1)*iKr(:,els).* ...
                jdh_(:,els).*YYp - kk*(kk+1)*iKr_(:,els).* ...
                djh_(:,els).*YpY)),1));

            if ~pa.Results.internal
              RgJ11(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(:,els).* ...
                  rM(:,els).*jj_(:,els).*(YtYp-YpYt),1));
              RgJ22(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(:,els).* ...
                  (rM(:,els).*djdj_(:,els).*(YtYp-YpYt) ...
                  -thetaM(:,els).*(jj*(jj+1)*iKr(:,els).* ...
                  jdj_(:,els).*YYp - kk*(kk+1)*iKr_(:,els).* ...
                  djj_(:,els).*YpY)),1));
            end
          end

          if ~(mirrorsym&&(rem(jj,2)~=rem(kk,2)))
            J12(jindx) = mfac.*(pi2*Nm(jj,kk)*sum(dS(:,els).* ...
                (rM(:,els).*jdh_(:,els).*(YpYp+YtYt) ...
                -kk*(kk+1)*iKr_(:,els).*thetaM(:,els).* ...
                jh_(:,els).*YtY),1));
            J21(jindx) = mfac.*(pi2*Nm(jj,kk)*sum(dS(:,els).* ...
                (-rM(:,els).*djh_(:,els).*(YtYt+YpYp) ...
                +jj*(jj+1)*iKr(:,els).*thetaM(:,els).*jh_(:,els).*YYt),1));

            if ~pa.Results.internal
              RgJ12(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(:,els).* ...
                  (rM(:,els).*jdj_(:,els).*(YpYp+YtYt) ...
                  -kk*(kk+1)*iKr_(:,els).*thetaM(:,els).*jj_(:,els).*YtY),1));
              RgJ21(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(:,els).* ...
                  (-rM(:,els).*djj_(:,els).*(YtYt+YpYp) ...
                  +jj*(jj+1)*iKr(:,els).*thetaM(:,els).*jj_(:,els).*YYt),1));
            end
          end
          jcount=jcount+length(indexzj);

        end
      end

      %Set up the Q and RgQ for sparse input:
      Qv(1:lengthJs)=1i*k_medium*(k_particle*J21+k_medium*J12);
      Qv(lengthJs+1:lengthJs*2)=1i*k_medium*(k_particle*J11+k_medium*J22);
      Qv(2*lengthJs+1:lengthJs*3)=1i*k_medium*(k_particle*J22+k_medium*J11);
      Qv(3*lengthJs+1:lengthJs*4)=1i*k_medium*(k_particle*J12+k_medium*J21);

      if ~pa.Results.internal
        RgQv(1:lengthJs)=1i*k_medium*(k_particle*RgJ21+k_medium*RgJ12);
        RgQv(lengthJs+1:lengthJs*2)=1i*k_medium* ...
            (k_particle*RgJ11+k_medium*RgJ22);
        RgQv(2*lengthJs+1:lengthJs*3)=1i*k_medium* ...
            (k_particle*RgJ22+k_medium*RgJ11);
        RgQv(3*lengthJs+1:lengthJs*4)=1i*k_medium* ...
            (k_particle*RgJ12+k_medium*RgJ21);
      end

      %Convert the vectors of Q and RgQ into sparse vector form:
      i_s=[i_indexes;i_indexes;i_indexes+Nmax*(Nmax+2); ...
          i_indexes+Nmax*(Nmax+2)];
      j_s=[j_indexes;j_indexes+Nmax*(Nmax+2);j_indexes; ...
          j_indexes+Nmax*(Nmax+2)];
        
      Q=sparse(i_s,j_s,Qv,2*Nmax*(Nmax+2),2*Nmax*(Nmax+2));
      if pa.Results.verbose
        disp(['Q cond: ' num2str(cond(Q))]);
      end
      if ~all(isfinite(Q(:)))
        error('Q is not finite');
      end
      
      if ~pa.Results.internal
        RgQ=sparse(i_s,j_s,RgQv,2*Nmax*(Nmax+2),2*Nmax*(Nmax+2));
        if pa.Results.verbose
          disp(['RgQ cond: ' num2str(cond(RgQ))]);
        end
        if ~all(isfinite(RgQ(:)))
          error('RgQ is not finite');
        end
      end

      if nargout > 1
        internal = tmatrix.setType('internal');
        internal.data = -inv(Q);
      end

        %solve the T-matrix:
        switch invmethod
          case 'forwardslash'
            % If we trust matlab or may be good for conductive particles
            tmatrix.data=-(RgQ/Q);
          case 'backslash'
            % If we trust matlab or may be good for conductive particles
            tmatrix.data=-(Q'\RgQ')';
          case 'inv'
            % This method probably uses LU decomposition, works better for
            % non-conductive and weakly conductive particles.
            % Wielaard, et al., Appl. Opt. 36, 4305-4313 (1997)
            tmatrix.data=-(RgQ * inv(Q));
          case 'pinv'
            tmatrix.data=-(RgQ * pinv(full(Q)));
          otherwise
            error('Unknown matrix inversion method');
        end

        % Store the type of T-matrix
        tmatrix.setType('scattered');
      end
    end
  end

  methods (Hidden)
    function shape = getGeometry(tmatrix, wavelength)
      % Get AxisymInterp shape with radius in specified units

      % Convert points to cylindrical coordinates
      rhoz = tmatrix.points(1, :).* ...
          [sin(tmatrix.points(2, :)); cos(tmatrix.points(2, :))];

      shape = ott.shapes.AxisymInterp(rhoz.*wavelength, ...
          'position', tmatrix.position*wavelength, ...
          'rotation', tmatrix.rotation);
    end
  end

  methods % Getters/setter
    function tmatrix = set.points(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 2, ...
          'points must be 2xN numeric matrix');
      assert(all(val(2, :) >= 0 & val(2, :) <= pi), ...
          'points(2, :) must all be in range [0, pi]');
      assert(all(val(1, :) >= 0), ...
          'points(1, :) must be positive');
      tmatrix.points = val;
    end

    function tmatrix = set.normals(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 2, ...
          'normals must be 2xN numeric matrix');
      tmatrix.normals = val;
    end

    function tmatrix = set.areas(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 1, ...
          'areas must be 1xN numeric matrix');
      tmatrix.areas = val;
    end
  end
end
