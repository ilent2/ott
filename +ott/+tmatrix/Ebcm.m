classdef Ebcm < ott.tmatrix.Tmatrix
% Constructs a T-matrix using extended boundary conditions method.
% Inherits from :class:`ott.tmatrix.Tmatrix`.
%
% Implements the extended boundary conditions methods for rotationally
% symmetric homogeneous particles.
%
% Properties
%   - points         -- (2xN numeric) Surface points [r; theta]
%   - normals        -- (2xN numeric) Normals at points [nr; ntheta]
%   - areas          -- (1xN numeric) Conic section surface areas
%   - relative_index -- (numeric) Relative refractive index of particle.
%   - xySymmetry     -- (logical) True if using xy-mirror optimisations
%   - invMethod      -- Inversion method used for T-matrix calculation
%
% Static methods
%   - FromShape   -- Construct from geometric shape object.
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
    relative_index % (numeric) Relative refractive index of particle
    xySymmetry     % (logical) True if using xy-mirror optimisations
    invMethod      % Inversion method used for T-matrix calculation
  end

  methods (Static)
    function varargout = FromShape(varargin)
      % Construct a T-matrix using EBCM from a shape object.
      %
      % If the shape object is not of type :class:`ott.shape.AxisymLerp`,
      % first casts the shape to this type using the default cast (if
      % supported, otherwise raises an error).
      %
      % Usage
      %   tmatrix = Ebcm.FromAxisymInterpShape(shape, relative_index, ...)
      %   Calculate external T-matrix (unless `internal` is true)
      %
      %   [external, internal] = Ebcm.FromAxisymInterpShape(...)
      %   Calculate both the internal and external T-matrices.
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- Description of shape geometry.
      %     Object must be a AxisymInterp or be castable to AxisymInterp.
      %
      %   - relative_index (numeric) -- Particle relative refractive index.
      %
      % See :meth:`Ebcm` for additional named parameters.

      p = inputParser;
      p.addRequired('shape', @(x) isa(x, 'ott.shape.Shape'));
      p.addOptional('relative_index', 1.0);
      p.addParameter('Nmax', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Ensure shape is an AxisymInterp
      shape = p.Results.shape;
      xySymmetry = shape.xySymmetry;
      if ~isa(shape, 'ott.shape.AxisymInterp')
        shape = ott.shape.AxisymInterp(shape);
      end

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
      [points, normals, ds] = shape.boundaryPoints(npts);

      % Construct T-matrices
      [varargout{1:nargout}] = ott.tmatrix.Ebcm(points, normals, ds, ...
          'relative_index', p.Results.relative_index, 'Nmax', Nmax, ...
          'xySymmetry', xySymmetry, unmatched{:});
    end
  end

  methods
    function [Texternal, Tinternal, data] = Ebcm(varargin)
      % Calculates T-matrix using extended boundary condition method
      %
      % Usage
      %   tmatrix = Ebcm(points, normals, area, relative_index, ...)
      %   Calculate external T-matrix (unless `internal` is true)
      %
      %   [external, internal, data] = Ebcm(...)
      %   Calculate both the internal and external T-matrices, and return
      %   the VswfData structure used during calculations.
      %
      % Parameters
      %   - points (2xN numeric) -- Coordinates describing surface.
      %     Spherical coordinates (omitting azimuthal angle: `[r; theta]`).
      %
      %   - normals (2xN numeric) -- Normals at surface points.
      %
      %   - areas (N numeric) -- Area elements at surface points.
      %
      %   - relative_index (numeric) -- Particle relative refractive index.
      %
      % Optional named parameters
      %   - xySymmetry (logical) -- If calculation should use xy-mirror
      %     symmetry optimisations.  Default: ``false``.
      %     Doesn't check points describe valid mirror symmetric shape.
      %
      %   - invMethod (enum|function handle) -- Inversion method for
      %     T-matrix calculation.  Currently supported methods are
      %     'forwardslash', 'backslash', 'inv', 'pinv'.  A custom
      %     inversion method can be specified using a function handle
      %     with the format ``inv_method(RgQ, Q)`` which returns the
      %     T-matrix data.  Default: ``forwardslash``.  Only used for
      %     external T-matrix calculation, may change in future.
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
      %     of the Q and RgQ matrices.  Default: ``false``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      p = inputParser;
      p.addRequired('points', @isnumeric);
      p.addRequired('normals', @isnumeric);
      p.addRequired('areas', @isnumeric);
      p.addOptional('relative_index', 1.0, @isnumeric);
      p.addParameter('xySymmetry', false);
      p.addParameter('invMethod', 'forwardslash');
      p.addParameter('internal', false);
      p.addParameter('Nmax', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.addParameter('verbose', []);
      p.parse(varargin{:});

      % Store parameters (changing units)
      Texternal.points = p.Results.points;
      Texternal.normals = p.Results.normals;
      Texternal.areas = p.Results.areas;
      Texternal.relative_index = p.Results.relative_index;
      Texternal.xySymmetry = p.Results.xySymmetry;
      Texternal.invMethod = p.Results.invMethod;

      calcExternalData = nargout == 2 || (nargout == 1 && ~p.Results.internal);

      % Get or estimate Nmax from the inputs
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        maxRadius = max(Texternal.points(1, :));
        if p.Results.internal
          Nmax = ott.utils.ka2nmax(maxRadius*2*pi*Texternal.index_relative);
        else
          Nmax = ott.utils.ka2nmax(maxRadius*2*pi);
        end
      else
        assert(isnumeric(Nmax) && isscalar(Nmax) && Nmax > 0, ...
            'Nmax must be positive numeric scalar');
      end

      % Get kr inside/outside particle
      kr=2*pi*Texternal.relative_index*Texternal.points(1, :);
      ikr=1./kr;
      kr_=2*pi*Texternal.points(1, :);
      ikr_=1./kr_;

      % Pre-calculate spherical harmonics and Bessel functions
      data = p.Results.data;
      ci = 1:ott.utils.combined_index(Nmax, Nmax);
      [Y, Ytheta, Yphi, data] = data.evaluateYtp(ci, ...
          Texternal.points(2, :), 0);
      [jkr, djkr, data] = data.evaluateBessel(1:Nmax, kr, 'regular');
      [jkr_, djkr_, data] = data.evaluateBessel(1:Nmax, kr_, 'regular');
      [hkr_, dhkr_, data] = data.evaluateBessel(1:Nmax, kr_, 'outgoing');
      hkr_ = hkr_ * 2;
      dhkr_ = dhkr_ * 2;

      Nn=1./sqrt((1:Nmax).*((1:Nmax)+1));
      Nm = Nn.' .* Nn;

      % Allocate memory for Q/RgQ matrix data
      lengthJs=Nmax^2*(Nmax+2)-sum((1:Nmax-1).*((1:Nmax-1)+1));
      J11=zeros(lengthJs,1);
      J12=J11;
      J21=J11;
      J22=J11;
      if calcExternalData
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

      % TODO: Do we still need these repmats?  Or is it better to just
      % use bsxfun style operations?
      fillnum=(2*Nmax+1);
      rM=repmat(Texternal.normals(1,:),[fillnum,1]);
      thetaM=repmat(Texternal.normals(2,:),[fillnum,1]);
      iKr=repmat(ikr,[fillnum,1]);
      iKr_=repmat(ikr_,[fillnum,1]);
      dS=repmat(Texternal.areas,[fillnum,1]);

      pi2=2*pi;

      k_medium = 2*pi;
      k_particle = 2*pi*Texternal.relative_index;

      for jj=1:Nmax
        %remember that when phi = 0, conj(Yphi)=-Yphi, conj(Ytheta)=Ytheta.

        for kk=1:Nmax

          px=min(kk,jj);
          pxrng = -px:px;

          els= 1:length(pxrng);

          jindx=jcount+els;

          i_indexes(jindx)=kk*(kk+1)+pxrng;
          j_indexes(jindx)=jj*(jj+1)+pxrng;

          jci = ott.utils.combined_index(jj, pxrng);
          kci = ott.utils.combined_index(kk, pxrng);

          %Calculate the cross products of the spherical harmonics
          YtYp=-Ytheta(jci, :).*Yphi(kci, :);
          YpYt=Yphi(jci, :).*Ytheta(kci, :);

          YYt=Y(jci, :).*Ytheta(kci, :);
          YtY=Ytheta(jci, :).*Y(kci, :);

          YYp=-Y(jci, :).*Yphi(kci, :);
          YpY=Yphi(jci, :).*Y(kci, :);

          YpYp=-Yphi(jci, :).*Yphi(kci, :);
          YtYt=Ytheta(jci, :).*Ytheta(kci, :);

          %calculate the cross products of spherical bessel functions.
          jh_=jkr(jj, :).*hkr_(kk, :);
          if calcExternalData
            jj_=jkr(jj, :).*jkr_(kk, :);
          end

          jdh_=jkr(jj, :).*dhkr_(kk, :);
          djh_=djkr(jj, :).*hkr_(kk, :);

          if calcExternalData
            jdj_=jkr(jj, :).*djkr_(kk, :);
            djj_=djkr(jj, :).*jkr_(kk, :);
          end

          djdh_=djkr(jj, :).*dhkr_(kk, :);
          if calcExternalData
            djdj_=djkr(jj, :).*djkr_(kk, :);
          end

          % perform the cross product followed by dotting the
          % normal vector and summing over the area elements.

          mfac=1;
          if ~(Texternal.xySymmetry && (rem(jj,2)==rem(kk,2)))
            J11(jindx) = mfac.*(pi2*Nm(jj,kk)* ...
                sum(dS(els,:).*rM(els,:).*jh_.*(YtYp-YpYt),2));
            J22(jindx) = mfac.*(pi2*Nm(jj,kk)* ...
                sum(dS(els,:).*(rM(els,:).*djdh_.*(YtYp-YpYt) ...
                -thetaM(els,:).*(jj*(jj+1)*iKr(els,:).* ...
                jdh_.*YYp - kk*(kk+1)*iKr_(els,:).* ...
                djh_.*YpY)),2));

            if calcExternalData
              RgJ11(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(els,:).* ...
                  rM(els,:).*jj_.*(YtYp-YpYt),2));
              RgJ22(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(els,:).* ...
                  (rM(els,:).*djdj_.*(YtYp-YpYt) ...
                  -thetaM(els,:).*(jj*(jj+1)*iKr(els,:).* ...
                  jdj_.*YYp - kk*(kk+1)*iKr_(els,:).* ...
                  djj_.*YpY)),2));
            end
          end

          if ~(Texternal.xySymmetry && (rem(jj,2)~=rem(kk,2)))
            J12(jindx) = mfac.*(pi2*Nm(jj,kk)*sum(dS(els,:).* ...
                (rM(els,:).*jdh_.*(YpYp+YtYt) ...
                -kk*(kk+1)*iKr_(els,:).*thetaM(els,:).* ...
                jh_.*YtY),2));
            J21(jindx) = mfac.*(pi2*Nm(jj,kk)*sum(dS(els,:).* ...
                (-rM(els,:).*djh_.*(YtYt+YpYp) ...
                +jj*(jj+1)*iKr(els,:).*thetaM(els,:).*jh_.*YYt),2));

            if calcExternalData
              RgJ12(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(els,:).* ...
                  (rM(els,:).*jdj_.*(YpYp+YtYt) ...
                  -kk*(kk+1)*iKr_(els,:).*thetaM(els,:).*jj_.*YtY),2));
              RgJ21(jindx)=mfac.*(pi2*Nm(jj,kk)*sum(dS(els,:).* ...
                  (-rM(els,:).*djj_.*(YtYt+YpYp) ...
                  +jj*(jj+1)*iKr(els,:).*thetaM(els,:).*jj_.*YYt),2));
            end
          end
          jcount=jcount+length(pxrng);

        end
      end

      % Get indices for sparse Q/RgQ packing
      i_s=[i_indexes;i_indexes;i_indexes+Nmax*(Nmax+2); ...
          i_indexes+Nmax*(Nmax+2)];
      j_s=[j_indexes;j_indexes+Nmax*(Nmax+2);j_indexes; ...
          j_indexes+Nmax*(Nmax+2)];

      % Package Q for inversion
      Qv(1:lengthJs)=1i*k_medium*(k_particle*J21+k_medium*J12);
      Qv(lengthJs+1:lengthJs*2)=1i*k_medium*(k_particle*J11+k_medium*J22);
      Qv(2*lengthJs+1:lengthJs*3)=1i*k_medium*(k_particle*J22+k_medium*J11);
      Qv(3*lengthJs+1:lengthJs*4)=1i*k_medium*(k_particle*J12+k_medium*J21);

      Q=sparse(i_s,j_s,Qv,2*Nmax*(Nmax+2),2*Nmax*(Nmax+2));
      if p.Results.verbose
        disp(['Q cond: ' num2str(cond(Q))]);
      end
      if ~all(isfinite(Q(:)))
        error('Q is not finite');
      end

      % Package RgQ for inversion
      if calcExternalData
        RgQv(1:lengthJs)=1i*k_medium*(k_particle*RgJ21+k_medium*RgJ12);
        RgQv(lengthJs+1:lengthJs*2)=1i*k_medium* ...
            (k_particle*RgJ11+k_medium*RgJ22);
        RgQv(2*lengthJs+1:lengthJs*3)=1i*k_medium* ...
            (k_particle*RgJ22+k_medium*RgJ11);
        RgQv(3*lengthJs+1:lengthJs*4)=1i*k_medium* ...
            (k_particle*RgJ12+k_medium*RgJ21);

        RgQ=sparse(i_s,j_s,RgQv,2*Nmax*(Nmax+2),2*Nmax*(Nmax+2));
        if p.Results.verbose
          disp(['RgQ cond: ' num2str(cond(RgQ))]);
        end
        if ~all(isfinite(RgQ(:)))
          error('RgQ is not finite');
        end
      end

      % Calculate external T-matrix
      if calcExternalData

        if isa(Texternal.invMethod, 'function_handle')

          Texternal.data = Texternal.invMethod(RgQ, Q);

        else

          switch Texternal.invMethod
            case 'forwardslash'
              % If we trust matlab or may be good for conductive particles
              Texternal.data=-(RgQ/Q);
            case 'backslash'
              % If we trust matlab or may be good for conductive particles
              Texternal.data=-(Q'\RgQ')';
            case 'inv'
              % This method probably uses LU decomposition, works better for
              % non-conductive and weakly conductive particles.
              % Wielaard, et al., Appl. Opt. 36, 4305-4313 (1997)
              Texternal.data=-(RgQ * inv(Q)); %#ok<MINV>
            case 'pinv'
              Texternal.data=-(RgQ * pinv(full(Q)));
            otherwise
              error('Unknown matrix inversion method');
          end
        end

        Texternal = Texternal.setType('scattered');
      end

      % Calculate internal T-matrix from Q
      if (p.Results.internal && nargout == 1) || nargout == 2

        Tinternal = Texternal;
        Tinternal.data = -inv(Q);
        Tinternal = Tinternal.setType('internal');
        Tinternal.invMethod = 'inv';

        if nargout == 1
          Texternal = Tinernal;
        end
      end
    end
  end

  methods % Getters/setter
    function tmatrix = set.relative_index(tmatrix, val)
      assert(isnumeric(val) && isscalar(val), ...
          'relative index must be numeric scalar');
      tmatrix.relative_index = val;
    end

    function tmatrix = set.points(tmatrix, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 2, ...
          'points must be 2xN numeric matrix');
      assert(all(val(2, :) >= 0 & val(2, :) <= pi), ...
          'points(2, :) must all be in range [0, pi]');
      assert(all(val(1, :) > 0), ...
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
