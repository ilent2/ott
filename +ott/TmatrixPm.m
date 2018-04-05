classdef TmatrixPm < ott.Tmatrix
%TmatrixPm constructs a T-matrix using the point matching method
%
% TmatrixPm properties:
%   k_medium          Wavenumber in the surrounding medium
%   k_particle        Wavenumber of the particle
%
% TmatrixPm methods:
%   getInternal       Get the internal T-matrix
%
% This class is based on tmatrix_pm.m from ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    k_medium            % Wavenumber of medium
    k_particle          % Wavenumber of particle

    idata               % Internal T-matrix
  end

  methods (Static)
    function tmatrix = simple(shape, parameters, varargin)
      %SIMPLE construct a T-matrix using PM for a simple shape.
      %
      % Supported shapes [parameters]:
      %   'ellipsoid'       Ellipsoid [ a b c]
      %   'cylinder'        z-axis aligned cylinder [ radius height ]
      %   'superellipsoid'  Superellipsoid [ a b c e n ]
      %   'cone-tipped-cylinder'      [ radius height cone_height ]
      %   'cube'            Cube [ width ]
      %
      %  TMATRIXPM(..., 'Nmax', Nmax) specifies the size of the
      %  T-matrix to use.  If not specified, the size is calculated
      %  from ott.utils.ka2nmax(max_radius*k_medium).
      %
      %  TMATRIXPM(..., 'k_medium', k)
      %  or TMATRIXPM(..., 'wavelength_medium', wavelength)
      %  or TMATRIXPM(..., 'index_medium', index)
      %  specify the wavenumber, wavelength or index in the medium.
      %
      %  TMATRIXPM(..., 'wavelength0', wavelength) specifies the
      %  wavelength in the vecuum, required when index_particle or
      %  index_medium are specified.

      % Handle different shapes
      if strcmp(shape, 'ellipsoid')
        shape_idx = 0;
        r_max = max(parameters);
      elseif strcmp(shape, 'cylinder')
        shape_idx = 1;
        r_max = sqrt(parameters(1)^2 + parameters(2)^2/4.0);
      elseif strcmp(shape, 'superellipsoid')
        shape_idx = 2;
        r_max = max(parameters(1:3));
      elseif strcmp(shape, 'cone-tipped-cylinder')
        shape_idx = 3;
        r_max = sqrt(parameters(1)^2 + parameters(2)^2/4.0);
        r_max = max(r_max, parameters(2)/2.0 + parameters(3));
      elseif strcmp(shape, 'cube')
        shape_idx = 4;
        r_max = sqrt(3.0*parameters(1)^2/4.0);
      else
        error('Unsupported shape type');
      end

      % Parse optional parameters
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('Nmax', []);
      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength0', []);
      p.parse(varargin{:});

      % Get or estimate Nmax from the inputs
      if isempty(p.Results.Nmax)
        k_medium = ott.Tmatrix.parser_k_medium(p);
        Nmax = ka2nmax(r_max * k_medium);
      else
        Nmax = p.Results.Nmax;
      end

      % Determine if we have rotational symetry
      [~,~,rotational_symmetry] = shapesurface([],[],shape_idx,parameters);

      if rotational_symmetry
        ntheta = 4*(Nmax + 2);
        nphi = 1;
      else
        ntheta = 2*(Nmax + 2);
        nphi = 3*(Nmax + 2)+1;
      end

      [theta,phi] = angulargrid(ntheta,nphi);

      % TODO: Different ways to distribute points (random points)

      [r,normals] = shapesurface(theta,phi,shape_idx,parameters);

      % TODO: What does inputParser do with duplicate parameters?
      ott.TmatrixPm(r, theta, phi, normals, varargin{:}, 'Nmax', Nmax, ...
          'rotational_symmetry', rotational_symmetry);
    end
  end

  methods
    function tmatrix = TmatrixPm(r, theta, phi, normals, varargin)
      %TMATRIXPM calculates T-matrix using the point matching method
      %
      % TMATRIXPM(r, theta, phi, normals) uses points at r, theta, phi
      % with normals to calculate the T-matrix.
      %     r     radius of the point
      %     theta polar angle from +z axis (rad)
      %     phi   azimuthal angle, measured from +x towards +y axes (rad)
      %
      % Both the external and internal T-matricies are calculated.
      % The external T-matrix is encapsulated by this object, the
      % internal T-matrix can be retrieved using getInternal.
      % If 'internal' is requested using the optional parameter (see bellow),
      % the external T-matrix is replaced with the internal T-matrix.
      %
      %  TMATRIXPM(..., 'Nmax', Nmax) specifies the size of the
      %  T-matrix to use.  If not specified, the size is calculated
      %  from ott.utils.ka2nmax(max_radius*k_medium).
      %
      %  TMATRIXPM(..., 'k_medium', k)
      %  or TMATRIXPM(..., 'wavelength_medium', wavelength)
      %  or TMATRIXPM(..., 'index_medium', index)
      %  specify the wavenumber, wavelength or index in the medium.
      %
      %  TMATRIXPM(..., 'k_particle', k)
      %  or TMATRIXPM(..., 'wavelength_particle', wavelength)
      %  or TMATRIXPM(..., 'index_particle', index)
      %  specify the wavenumber, wavelength or index in the particle.
      %
      %  TMATRIXPM(..., 'wavelength0', wavelength) specifies the
      %  wavelength in the vecuum, required when index_particle or
      %  index_medium are specified.
      %
      %  TMATRIXPM(..., 'internal', internal) if true, encapsulates
      %  the internal T-matrix and discards the external T-matrix.
      %
      %  TMATRIXPM(..., 'rotational_symmetry', sym) if true the particle
      %  is assumed to be rotationally symmetric about the z-axis and
      %  phi must be all the same angle.

      tmatrix = tmatrix@ott.Tmatrix();

      % Parse inputs
      p = inputParser;
      p.addParameter('Nmax', []);
      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('k_particle', []);
      p.addParameter('wavelength_particle', []);
      p.addParameter('index_particle', []);
      p.addParameter('wavelength0', []);
      p.addParameter('internal', false);
      p.addParameter('rotational_symmetry', false);
      p.parse(varargin{:});

      % Store inputs k_medium and k_particle
      tmatrix.k_medium = tmatrix.parser_k_medium(p);
      tmatrix.k_particle = tmatrix.parser_k_particle(p);

      % Get or estimate Nmax from the inputs
      if isempty(p.Results.Nmax)
        Nmax = ka2nmax(max(r) * tmatrix.k_medium);
      else
        Nmax = p.Results.Nmax;
      end

      npoints = length(theta);
      total_orders = Nmax * (Nmax+2);

      % 3 vector components at each point, c/d,p/q coefficient per order
      coeff_matrix = zeros(6*npoints,4*total_orders);
      incident_wave_matrix = zeros(6*npoints,2*total_orders);

      if p.Results.rotational_symmetry
         T = sparse(2*total_orders,2*total_orders);
         T2 = sparse(2*total_orders,2*total_orders);
      else
         T = zeros(2*total_orders,2*total_orders);
         T2 = zeros(2*total_orders,2*total_orders);
      end

      import ott.utils.*

      k_relative = k_particle/k_medium;

      for n = 1:Nmax
        for m = -n:n

          % INCIDENT-SCATTERED
          [M1,N1,~,~,M2,N2] = vswf(n,m,tmatrix.k_medium*r,theta,phi);
          [M3,N3] = vswf(n,m,tmatrix.k_particle*r,theta,phi,3);

          ci = combined_index(n,m);

          M1 = perpcomponent(M1,normals);
          N1 = perpcomponent(N1,normals);
          M2 = perpcomponent(M2,normals);
          N2 = perpcomponent(N2,normals);
          M3 = perpcomponent(M3,normals);
          N3 = perpcomponent(N3,normals);
          M1 = M1(:);
          N1 = N1(:);
          M2 = M2(:);
          N2 = N2(:);
          M3 = M3(:);
          N3 = N3(:);

          % 1 is outgoing field, 3 is particle field, 2 is incoming field
          coeff_matrix(:,ci) = - [ M1; N1 ];
          coeff_matrix(:,ci+total_orders) = - [ N1; M1 ];
          coeff_matrix(:,ci+2*total_orders) = [ M3; k_relative*N3 ];
          coeff_matrix(:,ci+3*total_orders) = [ N3; k_relative*M3 ];

          incident_wave_matrix(:,ci) = [ M2; N2 ];
          incident_wave_matrix(:,ci+total_orders) = [ N2; M2 ];

          % TODO: Output progress
        end
      end

      for n = 1:Nmax
        for m = -n:n

          ci = combined_index(n,m);

          if p.Results.rotational_symmetry
            number_of_nm = 1 + Nmax - max(abs(m),1);
            nm_to_use = combined_index(max(abs(m),1):Nmax, ...
                ones(1,number_of_nm)*m);
            nm_to_use = [ nm_to_use nm_to_use+total_orders ];
            all_indices = [ nm_to_use  nm_to_use+2*total_orders ];

            incident_wave_vector = incident_wave_matrix(:,ci);
            Tcol = coeff_matrix(:,all_indices) \ incident_wave_vector;
            T(nm_to_use,ci) = Tcol(1:(2*number_of_nm),1);
            T2(nm_to_use,ci) = Tcol((1+2*number_of_nm):(4*number_of_nm),1);

            incident_wave_vector = incident_wave_matrix(:,ci+total_orders);
            Tcol = coeff_matrix(:,all_indices) \ incident_wave_vector;
            T(nm_to_use,ci+total_orders) = Tcol(1:(2*number_of_nm),1);
            T2(nm_to_use,ci+total_orders) = ...
                Tcol((1+2*number_of_nm):(4*number_of_nm),1);
          else
            incident_wave_vector = incident_wave_matrix(:,ci);
            Tcol = coeff_matrix \ incident_wave_vector;
            T(:,ci) = Tcol(1:2*total_orders,1);
            T2(:,ci) = Tcol((1+2*total_orders):4*total_orders,1);

            incident_wave_vector = incident_wave_matrix(:,ci+total_orders);
            Tcol = coeff_matrix \ incident_wave_vector;
            T(:,ci+total_orders) = Tcol(1:2*total_orders,1);
            T2(:,ci+total_orders) = Tcol((1+2*total_orders):4*total_orders,1);
          end

          % TODO: Output progress
        end
      end

      % Store the T-matrix data
      if p.Results.internal
        tmatrix.data = T2;
        tmatrix.idata = [];
      else
        tmatrix.data = T;
        tmatrix.idata = T2;
      end

      % TODO: Store the type of T-matrix (scattered/total)
    end

    function tmatrix2 = getInternal(tmatrix)
      %GETINTERNAL get a T-matrix object for the internal T-matrix data

      if isempty(tmatrix.idata)
        error('Internal T-matrix data has been cleaned');
      end

      % TODO: We should probably store that the matrix is internal
      tmatrix2 = ott.Tmatrix(tmatrix.idata);
    end

    function cleanInternal(tmatrix)
      %CLEANINTERNAL remove the internal T-matrix data
      tmatrix.idata = [];
    end
  end
end
