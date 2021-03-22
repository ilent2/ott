classdef TmatrixMie < ott.Tmatrix
%TmatrixMie construct T-matrix from Mie scattering coefficients
%
% TmatrixMie properties:
%   radius            The radius of the sphere the T-matrix represents
%   k_medium          Wavenumber in the trapping medium
%   k_particle        Wavenumber of the particle
%
% This class is based on tmatrix_mie.m and tmatrix_mie_layered.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    radius            % Radius of particle
    k_medium          % Wavenumber of medium
    k_particle        % Wavenumber of particle
    
    mu_relative       % Relative permeability
  end

  methods (Access=protected)
    function T = tmatrix_mie(tmatrix, Nmax, internal)
      %TMATRIX_MIE code from tmatrix_mie.m

      n=[1:Nmax];

      m = tmatrix.k_particle/tmatrix.k_medium;
      mu = tmatrix.mu_relative;

      r0 = tmatrix.k_medium * tmatrix.radius;
      r1 = tmatrix.k_particle * tmatrix.radius;

      indexing=ott.utils.combined_index(1:Nmax^2+2*Nmax)';

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      j0 = (sbesselj(n,r0)).';
      j1 = (sbesselj(n,r1)).';
      h0 = (sbesselh1(n,r0)).';
      j0d = (sbesselj(n-1,r0) - n.*sbesselj(n,r0)/r0).';
      j1d = (sbesselj(n-1,r1) - n.*sbesselj(n,r1)/r1).';
      h0d = (sbesselh1(n-1,r0) - n.*sbesselh1(n,r0)/r0).';

      if internal == false
        % These differ from some definitions of a,b (a = -b, b = -a)
        b = -( mu*j1d.*j0 - m*j0d.*j1 ) ./ ( mu*j1d.*h0 - m*h0d.*j1 );
        a = -( mu*j0d.*j1 - m*j1d.*j0 ) ./ ( mu*h0d.*j1 - m*j1d.*h0 );

        T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
            [a(indexing);b(indexing)]);
      else
        % Calculate internal T-matrix
        d = ( h0d.*j0 - j0d.*h0 ) ./ ( m*j1.*h0d - j1d.*h0 );
        c = ( j0d.*h0 - h0d.*j0 ) ./ ( m*j1d.*h0 - h0d.*j1 );
        T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
            [c(indexing);d(indexing)]);
      end
    end

    function T = tmatrix_mie_layered(tmatrix, Nmax, internal)
      %TMATRIX_MIE code from tmatrix_mie_layered.m

      k_layer=[tmatrix.k_particle,tmatrix.k_medium]; % Medium on outside
      radius=[tmatrix.radius,tmatrix.radius(end)]; % medium same as final layer
      n_layer=k_layer/2/pi;

      n = 1:Nmax; %n for nmax

      lastElement=length(k_layer);

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      %generate all special functions first:
      if length(tmatrix.k_particle)>1
          [jN,jNd] = sbesselj(n,[k_layer.*radius, ...
              k_layer(2:end).*radius(1:end-1)]);
          [hN,hNd] = sbesselh1(n,[k_layer.*radius, ...
              k_layer(2:end).*radius(1:end-1)]);
      else
          [jN,jNd] = sbesselj(n,k_layer(:).*radius(:));
          [hN,hNd] = sbesselh1(n,k_layer(:).*radius(:));
      end
      jN=jN.';
      hN=hN.';
      jNd=jNd.';
      hNd=hNd.';

      d1_N=jNd./jN;
      d3_N=hNd./hN;
      r_N=jN./hN;

      ha_0=1/n_layer(1);
      hb_0=1;

      for ii=1:length(tmatrix.k_particle)
          ha_n=ha_0;
          hb_n=hb_0;

          m_0 = n_layer(ii);

          d1_0=d1_N(:,ii);
          d3_0=d3_N(:,ii);
          r_0=r_N(:,ii);

          if ii>1
              m_n = n_layer(ii-1);

              d1_n=d1_N(:,lastElement+ii-1);
              d3_n=d3_N(:,lastElement+ii-1);
              r_n=r_N(:,lastElement+ii-1);
          else
              m_n=1;

              d1_n=0;
              d3_n=0;
              r_n=0;
          end

          g0=m_0*ha_n-m_n*d1_n;
          g1=m_0*ha_n-m_n*d3_n;
          g0b=m_n*hb_n-m_0*d1_n;
          g1b=m_n*hb_n-m_0*d3_n;
          q_0=r_n./r_0;
          ha_0=(g1.*d1_0-q_0.*g0.*d3_0)./(g1-q_0.*g0);
          hb_0=(g1b.*d1_0-q_0.*g0b.*d3_0)./(g1b-q_0.*g0b);

      end
      m_0 = n_layer(lastElement-1);
      m_1 = n_layer(lastElement);

      m=m_0/m_1;

      d1_1=d1_N(:,lastElement);
      d3_1=d3_N(:,lastElement);
      r_1=r_N(:,lastElement);

      % %TM/TE coeffs...
      al_1=r_1.*(ha_0-m*d1_1)./(ha_0-m*d3_1);
      bl_1=r_1.*(m*hb_0-d1_1)./(m*hb_0-d3_1);

      %swap the modes to TE/TM... and negatize.
      b = -al_1;
      a = -bl_1;

      %t-matrix indices...
      indexing=ott.utils.combined_index(1:Nmax^2+2*Nmax)';

      if internal == false

        T=sparse(1:2*(Nmax^2+2*Nmax),1:2*(Nmax^2+2*Nmax),[a(indexing);b(indexing)]);

      else

        r_0=(jN(:,lastElement)./jN(:,lastElement-1));

        d = r_0.*(d3_1 - d1_1 )  ./ (ha_0-m*d3_1);
        c = r_0.*(d3_1 - d1_1 )  ./ (m*hb_0 - d3_1);

        warning('ott:tmatrix_mie_layered:internalcoefficientwarning', ...
            ['The internal coefficients are for the outermost layer only...' ...
             ' the real ones are only defined for each layer.']);
        T=sparse(1:2*(Nmax^2+2*Nmax),1:2*(Nmax^2+2*Nmax),[c(indexing);d(indexing)]);
      end
    end
  end
  
  methods (Static, Hidden)
    
    function p = input_parser(varargin)
      % Helper for input parsing
      
      p = inputParser;

      p.addParameter('progress_callback', ...
        @ott.TmatrixMie.DefaultProgressCallback);
      p.addParameter('Nmax', []);
      p.addParameter('wavelength0', []);
      p.addParameter('internal', false);
      p.addParameter('shrink', true);

      p.addParameter('index_relative', []);
      p.addParameter('mu_relative', 1.0);

      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);

      p.addParameter('k_particle', []);
      p.addParameter('wavelength_particle', []);
      p.addParameter('index_particle', []);
      
      p.parse(varargin{:});
    end
    
    function DefaultProgressCallback()
      % Default progress callback function
      
      % TODO: Progress callback for this method
    end
    
  end

  methods (Static)
    function tmatrix = simple(shape, varargin)
      %SIMPLE construct a T-matrix using Mie solution for a sphere.
      %
      % SIMPLE(shape) constructs a new simple T-matrix for the given
      % ott.shapes.Shape object.  Shape must implement a maxRadius
      % method, the result of which is used as the sphere radius.
      %
      % SIMPLE(name, parameters) constructs a new T-matrix for the
      % shape described by the name and parameters.  The name must
      % be 'sphere' and the parameters must be the radius.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('parameters', []);
      p.parse(varargin{:});
      
      % Handle different shapes
      if isa(shape, 'ott.shapes.Shape') && isempty(p.Results.parameters)
        radius = shape.maxRadius;
      elseif ischar(shape) && strcmpi(shape, 'sphere') ...
          && ~isempty(p.Results.parameters)
        radius = p.Results.parameters;
        varargin = varargin(2:end);
      else
        error('Must input either Shape object or string and parameters');
      end
      
      % Genreate an instance of this object
      tmatrix = ott.TmatrixMie(radius, varargin{:});
    end
  end

  methods
    function tmatrix = TmatrixMie(radius, varargin)
      %TMATRIXMIE construct a new Mie T-matrix for a sphere with size radius.
      %
      % tmatirx = TmatrixMie(radius, ...)
      %
      % If radius and k_particle (see bellow) are vectors, calculates the
      % coefficients for a leyered sphere with radius and k_particle
      % specifying the radius and wavenumber of each layer starting
      % from the core and going out.
      %
      % Note: Will not compute for particles of 500 layers or more. Not all of
      % the recursions in the reference are implemented because the original
      % author was lazy. For layered sphere, Nmax of 100 is numerically stable.
      %
      %  TMATRIXMIE(..., 'Nmax', Nmax) specifies the size of the
      %  T-matrix to use.  If not specified, the size is calculated
      %  from ott.utils.ka2nmax(radius*k_medium) or set to 100 for
      %  layered spheres.
      %
      %  TMATRIXMIE(..., 'k_medium', k)
      %  or TMATRIXMIE(..., 'wavelength_medium', wavelength)
      %  or TMATRIXMIE(..., 'index_medium', index)
      %  specify the wavenumber, wavelength or index in the medium.
      %  If not specified, the default is k_medium = 2*pi.
      %
      %  TMATRIXMIE(..., 'k_particle', k)
      %  or TMATRIXMIE(..., 'wavelength_particle', wavelength)
      %  or TMATRIXMIE(..., 'index_particle', index)
      %  specify the wavenumber, wavelength or index in the particle.
      %
      %  TMATRIXMIE(..., 'index_relative', n) can be used to specify
      %  the relative refractive index if either the medium or particle
      %  wavenumber can already be determined.
      %
      %  TMATRIXMIE(..., 'mu_relative', mu) specify the relative permeability.
      %  Caution: may change in future, doesn't adjust the refractive index.
      %
      %  TMATRIXMIE(..., 'wavelength0', wavelength) specifies the
      %  wavelength in the vecuum, required when index_particle or
      %  index_medium are specified.
      %
      %  TMATRIXMIE(..., 'internal', internal) if true, calculates the
      %  T-matrix for the internal coefficients.
      %
      % Layered sphere reference:
      % "Improved recursive algorithm for light scattering by a multilayered
      % sphere", Wen Yang, Applied Optics 42(9), 2003

      tmatrix = tmatrix@ott.Tmatrix();

      % Parse inputs
      p = ott.TmatrixMie.input_parser(varargin{:});

      % Store inputs: radius, k_medium, k_particle
      tmatrix.radius = radius;
      [tmatrix.k_medium, tmatrix.k_particle] = ...
          tmatrix.parser_wavenumber(p, 2.0*pi);
      tmatrix.mu_relative = p.Results.mu_relative;

      % Check radius and k_particle are similar lengths
      if numel(tmatrix.radius) ~= numel(tmatrix.k_particle)
        error('radius and k_particle must be the same length');
      end

      % Check number of layers
      if numel(tmatrix.radius) >= 100
        warning('May not work well for particles with >= 100 layers');
      end

      % If Nmax not specified, choose a good Nmax
      if isempty(p.Results.Nmax)

        % Different Nmax for internal and external
        if p.Results.internal == false
          Nmax = ott.utils.ka2nmax(tmatrix.k_medium*radius(end));
        else
          Nmax = ott.utils.ka2nmax(tmatrix.k_particle(end)*radius(end));
        end

      else
        Nmax = p.Results.Nmax;
      end

      % Store the T-matrix type
      if p.Results.internal == false
        tmatrix.type = 'scattered';
      else
        tmatrix.type = 'internal';
      end

      % Calculate the T-matrix and store it
      if numel(tmatrix.radius) == 1
        tmatrix.data = tmatrix.tmatrix_mie(Nmax, p.Results.internal);
      else

        % Layered method needs at least Nmax=100 (apparently)
        if p.Results.shrink && isempty(p.Results.Nmax)
          oldNmax = Nmax;
          Nmax = max(100, Nmax);
        end

        % Calculate T-matrix
        tmatrix.data = tmatrix.tmatrix_mie_layered(Nmax, p.Results.internal);

        % Shrink T-matrix (required type to be already set)
        if p.Results.shrink
          tmatrix.Nmax = oldNmax;
        end

      end
    end
  end
end
