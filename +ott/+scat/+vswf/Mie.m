classdef Mie < ott.scat.vswf.Tmatrix
% Construct T-matrix with Mie scattering coefficients
%
% The Mie coefficients describe the scattering of a sphere.
% They can also be used to give a reasonable estimate of the force
% for non-spherical particles when no other suitable method is available.
%
% This class supports calculating Mie coefficients for homogeneous
% spheres and layered spheres.
%
% The T-matrix coefficients are calculated when the constructor is
% called.  The class supports changing the properties after
% construction in order to set appropriate units, for example, to
% change the particle radius from wavelengths to meters
%
% .. code::matlab
%   tmatrix.wavelength_medium = 1.0e-6;
%
% This will result in the radius changing by a similar factor.
%
% ``permittivity`` and ``permeability`` are relative to the respective vacuum
% quantities.
%
% Layered sphere reference:
% "Improved recursive algorithm for light scattering by a multilayered
% sphere", Wen Yang, Applied Optics 42(9), 2003
%
% Properties
%   - radius         -- Radius of sphere or radii of sphere layers
%   - wavelength_medium -- Wavelength in the medium (default: 1.0)
%   - permittivity_relative  -- Relative permittivity (can be complex)
%   - permeability_relative  -- Relative permeability (default: 1.0)
%   - permittivity_medium -- Relative permittivity in medium (default: 1.0)
%   - permeability_medium -- Relative permeability in medium (default: 1.0)
%
% Dependent properties
%   - index_relative -- Relative refractive index of sphere or layers
%   - index_medium   -- Refractive index of surrounding medium
%   - index_particle -- Refractive index of sphere or layers
%   - wavenumber_medium   -- Wave-number of medium
%   - wavenumber_particle -- Wave-number of particle
%   - permittivity_particle -- Relative permittivity in particle
%   - permeability_particle -- Relative permeability in particle
%   - wavelength_particle -- Wavelength inside particle
%
% This class is based on tmatrix_mie.m and tmatrix_mie_layered.m from ottv1.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Hidden, SetAccess=protected)
    wavelengthInternal % Actual wavelength medium value
  end

  properties (SetAccess=protected)
    radius          % Radius of sphere or radii of sphere layers
    permittivity_relative   % Relative permittivity (can be complex)
    permeability_relative   % Relative permeability (default: 1.0)
  end

  properties
    permittivity_medium  % Permittivity in medium (default: 1.0)
    permeability_medium  % Permeability in medium (default: 1.0)
  end

  properties (Dependent)

    % These properties are stored internally
    % Changing these has no effect on the T-matrix data, it just
    % changes the length scales for radius/wavelength (co-dependent)
    wavelength_medium  % Wavelength in the medium (default: 1.0)

    % Entirely dependent properties
    index_relative  % Relative refractive index of sphere or layers
    index_medium    % Refractive index of surrounding medium
    index_particle  % Refractive index of sphere or layers
    wavenumber_medium    % Wave-number of medium
    wavenumber_particle  % Wave-number of particle
    permittivity_particle  % Permittivity in particle
    permeability_particle  % Permeability in particle
    wavelength_particle  % Wavelength inside particle
  end

  methods (Access=protected)
    function T = tmatrix_mie(tmatrix, Nmax, internal)
      %TMATRIX_MIE code from tmatrix_mie.m

      n=[1:Nmax];

      m = tmatrix.index_relative;
      mu = tmatrix.permeability_relative;

      r0 = tmatrix.wavenumber_medium * tmatrix.radius;
      r1 = tmatrix.wavenumber_particle * tmatrix.radius;

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

      % Medium on outside
      k_layer=[tmatrix.wavenumber_particle,tmatrix.wavenumber_medium];

      % medium same as final layer
      radius=[tmatrix.radius,tmatrix.radius(end)];

      n_layer=k_layer/2/pi;

      n = 1:Nmax; %n for nmax

      lastElement=length(k_layer);

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      %generate all special functions first:
      if length(tmatrix.wavenumber_particle)>1
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

      for ii=1:length(tmatrix.wavenumber_particle)
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

    function [wavelength_medium, permeability_relative, ...
        permittivity_relative, permittivity_medium, ...
        permeability_medium] = parseInputs(varargin)
      % Parse dependent inputs for class constructor

      p = ott.utils.RelatedArgumentParser;

      p.addRequired('wavelength_medium', 1.0);
      p.addRequired('permittivity_relative');
      p.addRequired('permeability_relative', 1.0);
      p.addRequired('permittivity_medium', 1.0);
      p.addRequired('permeability_medium', 1.0);

      p.addOptional('index_relative');
      p.addOptional('index_medium');
      p.addOptional('index_particle');
      p.addOptional('wavenumber_medium');
      p.addOptional('wavenumber_particle');
      p.addOptional('permittivity_particle');
      p.addOptional('permeability_particle');
      p.addOptional('wavelength_particle');

      % Add rules
      %p.addRule('index_relative = index_particle ./ index_medium');
      p.addRule('wavenumber_medium = 2*pi / wavelength_medium');
      p.addRule('wavenumber_particle = 2*pi / wavelength_particle');
      p.addRule(['index_relative = ', ...
          'sqrt(permittivity_relative * permeability_relative)']);
      p.addRule(['index_medium = ', ...
          'sqrt(permittivity_medium * permeability_medium)']);
      p.addRule(['index_particle = ', ...
          'sqrt(permittivity_particle * permeability_particle)']);
      p.addRule(['wavelength_medium .* index_medium = ', ...
          'wavelength_particle .* index_particle']);
      p.addRule(['permittivity_relative = ', ...
          'permittivity_particle ./ permittivity_medium']);
      p.addRule(['permeability_relative = ', ...
          'permeability_particle ./ permeability_medium']);

      % Parse inputs
      p.parse(varargin{:});

      % Assign outputs
      wavelength_medium = p.RequiredResults.wavelength_medium;
      permeability_relative = p.RequiredResults.permeability_relative;
      permittivity_relative = p.RequiredResults.permittivity_relative;
      permittivity_medium = p.RequiredResults.permittivity_medium;
      permeability_medium = p.RequiredResults.permeability_medium;
      
      assert(~isempty(permittivity_relative), ...
        'Unable to determine permittivity_relative from inputs');

    end

    function DefaultProgressCallback()
      % Default progress callback function

      % TODO: Progress callback for this method
    end
  end

  methods (Static)
    function tmatrix = simple(shape, varargin)
      % Construct a T-matrix using Mie solution for a sphere.
      %
      % Usage
      %   ott.scat.vswf.tmatrix.simple(shape)
      %   Constructs a new simple T-matrix for the given
      %   ott.shapes.Shape object.  Shape must implement a maxRadius
      %   method, the result of which is used as the sphere radius.
      %
      %   ott.scat.vswf.tmatrix.simple(name, parameters)
      %   Constructs a new T-matrix for the
      %   shape described by the name and parameters.  The name must
      %   be 'sphere' and the parameters must be the radius.
      %
      % Additional named arguments are passed to the Mie constructor.
      % See also Mie

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
      tmatrix = ott.scat.vswf.Mie(radius, varargin{:});
    end
  end

  methods
    function tmatrix = Mie(radius, varargin)
      % Construct a new Mie T-matrix for a sphere.
      %
      % To construct a Mie T-matrix the particle radius and relative
      % refractive index needs to be specified.  The relative refractive
      % index can be specified using any suitable combination of optional
      % arguments.
      %
      % If the radius rand refractive index are vectors, the class
      % calculates the Mie coefficients for a layered sphere.
      % Elements describe the sphere layers starting from the core.
      %
      % Note: Will not compute for particles of 500 layers or more. Not all of
      % the recursions in the reference are implemented because the original
      % author was lazy. For layered sphere, Nmax of 100 is numerically stable.
      %
      % Usage
      %   tmatirx = Mie(radius, ...)
      %
      % Parameters
      %   - radius (numeric) -- Radius of the sphere or radius of layers.
      %     Radius has the same units as the medium wavelength, i.e.
      %     ``radius / wavelength_medium`` is a dimensionless quantity.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Size of the VSWF expansion used
      %     for the T-matrix.  Can be reduced after construction.
      %     Default: ``ott.utis.ka2nmax(radius*wavenumber_medium)``
      %     or for layered sphere: ``100``.  If ``internal = true``,
      %     uses the particle wave-number for ``ka2nmax``.
      %
      %   - wavelength_medium -- Wavelength in the medium (default: 1.0)
      %   - permittivity_relative  -- Relative permittivity (can be complex)
      %   - permeability_relative  -- Relative permeability (default: 1.0)
      %   - permittivity_medium -- Permittivity in medium (default: 1.0)
      %   - permeability_medium -- Permeability in medium (default: 1.0)
      %
      %   - index_relative -- Relative refractive index of sphere or layers
      %   - index_medium   -- Refractive index of surrounding medium
      %   - index_particle -- Refractive index of sphere or layers
      %   - wavenumber_medium   -- Wave-number of medium
      %   - wavenumber_particle -- Wave-number of particle
      %   - permittivity_particle -- Permittivity in particle
      %   - permeability_particle -- Permeability in particle
      %   - wavelength_particle -- Wavelength inside particle
      %
      %   - internal (logical) -- If true, calculates the internal T-matrix.
      %     Default: ``false``.
      %   - shrink (logical) -- If true, attempts to reduce the T-matrix
      %     size after the calculation.  Default: ``numel(radius)>1``.
      %
      %   - progress_callback (function_handle) -- Function to call
      %     to report progress of the method.  Not yet implemented.

      tmatrix = tmatrix@ott.scat.vswf.Tmatrix();

      % Get radius from shape
      if isa(radius, 'ott.shapes.Shape')
        radius = radius.maxRadius;
      end

      % Store positional arguments
      tmatrix.radius = radius;

      % Parse optional arguments
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('Nmax', []);
      p.addParameter('internal', false);
      p.addParameter('shrink', numel(radius)>1);
      p.addParameter('progress_callback', ...
        @ott.scat.vswf.tmatrix.Mie.DefaultProgressCallback);
      p.parse(varargin{:});

      % Parse dependent arguments
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];
      [tmatrix.wavelengthInternal, tmatrix.permeability_relative, ...
          tmatrix.permittivity_relative, tmatrix.permittivity_medium, ...
          tmatrix.permeability_medium] = tmatrix.parseInputs(unmatched{:});

      % Get default Nmax
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        if numel(radius) > 1
          Nmax = 100;       % Recommended >= for stability
        else
          if p.Results.internal
            Nmax = ott.utils.ka2nmax(tmatrix.wavenumber_particle(end)*radius(end));
          else
            Nmax = ott.utils.ka2nmax(tmatrix.wavenumber_medium*radius(end));
          end
        end
      end

      % Store the T-matrix type
      if p.Results.internal
        tmatrix.type = 'internal';
      else
        tmatrix.type = 'scattered';
      end

      % Calculate the T-matrix and store it
      if numel(tmatrix.radius) == 1
        tmatrix.data = tmatrix.tmatrix_mie(Nmax, p.Results.internal);
      else
        tmatrix.data = tmatrix.tmatrix_mie_layered(Nmax, p.Results.internal);
      end
      
      % Apply shrink if requested (find modes with sigificant power)
      if p.Results.shrink
        
        % Threshold for shrinking
        tol = 1.0e-6;
        
        % Find small ab coefficients
        ab = reshape(diag(tmatrix.data), [], 2);
        ab_small = abs(ab) < tol;
        
        % Get index of last non-small coeff.
        not_small = ~(ab_small(:, 1) & ab_small(:, 2));
        I = find(not_small, 1, 'last');
        
        % Convert row index to Nmax
        [new_Nmax, ~] = ott.utils.combined_index(I);
        tmatrix.Nmax = min(new_Nmax + 1, tmatrix.Nmax);
      end
    end
  end

  methods

    function tmatrix = set.radius(tmatrix, val)

      % Allow the radius to change size on initial set
      assert(isnumeric(val) && size(val, 1) == 1, ...
        'radius must be numeric 1xN vector');

      % Check values are ascending
      assert(all(diff(val) > 0), 'Radius values must be ascending');

      tmatrix.radius = val;
    end

    function tmatrix = set.wavelength_medium(tmatrix, val)

      assert(isnumeric(val) && isscalar(val), ...
        'wavelength_medium must be numeric scalar');

      % Adjust radius
      ratio = val ./ tmatrix.wavelengthInternal;
      tmatrix.radius = tmatrix.radius .* ratio;

      % Store the new value
      tmatrix.wavelengthInternal = val;
    end
    function val = get.wavelength_medium(tmatrix)
      val = tmatrix.wavelengthInternal;
    end

    %
    % Property validation
    %

    function tmatrix = set.wavelengthInternal(tmatrix, val)

      % Validate again for internal use
      assert(isnumeric(val) && isscalar(val), ...
        'wavelength_medium must be numeric scalar');

      tmatrix.wavelengthInternal = val;
    end
    function tmatrix = set.permittivity_relative(tmatrix, val)
      assert(isnumeric(val) && (isscalar(val) ...
        || all(size(val) == size(tmatrix.radius))), ...
        'permittivity_relative must be numeric and scalar or same size as radius');
      tmatrix.permittivity_relative = val;
    end
    function tmatrix = set.permeability_relative(tmatrix, val)
      assert(isnumeric(val) && (isscalar(val) ...
        || all(size(val) == size(tmatrix.radius))), ...
        'permeability_relative must be numeric and scalar or same size as radius');
      tmatrix.permeability_relative = val;
    end
    function tmatrix = set.permittivity_medium(tmatrix, val)
      assert(isnumeric(val) && isscalar(val), ...
        'permittivity_medium must be numeric scalar');
      tmatrix.permittivity_medium = val;
    end
    function tmatrix = set.permeability_medium(tmatrix, val)
      assert(isnumeric(val) && isscalar(val), ...
        'permeability_medium must be numeric scalar');
      tmatrix.permeability_medium = val;
    end

    %
    % Dependent properties
    %

    function val = get.index_relative(tmatrix)
      val = sqrt(tmatrix.permittivity_relative ...
          .* tmatrix.permeability_relative);
    end
    function val = get.index_medium(tmatrix)
      val = sqrt(tmatrix.permittivity_medium ...
          .* tmatrix.permeability_medium);
    end
    function val = get.index_particle(tmatrix)
      val = sqrt(tmatrix.permittivity_particle ...
          .* tmatrix.permeability_particle);
    end
    function val = get.wavenumber_medium(tmatrix)
      val = 2*pi ./ tmatrix.wavelength_medium;
    end
    function val = get.wavenumber_particle(tmatrix)
      val = 2*pi ./ tmatrix.wavelength_particle;
    end
    function val = get.permittivity_particle(tmatrix)
      val = tmatrix.permittivity_relative .* tmatrix.permittivity_medium;
    end
    function val = get.permeability_particle(tmatrix)
      val = tmatrix.permeability_relative .* tmatrix.permeability_medium;
    end
    function val = get.wavelength_particle(tmatrix)
      val = tmatrix.wavelength_medium ./ tmatrix.index_relative;
    end
  end
end
