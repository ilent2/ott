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
      %   'sphere'          Sphere [ radius ]
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
      %
      %  TMATRIXPM(..., 'distribution', m) specifies point distriution method
      %      'angulargrid'    uses an angular grid of points
      %      'random'         uses randomly distributed points

      % Handle different shapes
      if ischar(shape)
        shape = ott.shapes.Shape.simple(shape, parameters);
      end

      % Check the particle is star shaped
      if ~isa(shape, 'ott.shapes.StarShape')
        error('Only star shaped particles supported for now');
      end

      % Parse optional parameters
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('Nmax', []);
      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength0', []);
      p.addParameter('distribution', 'angulargrid');
      p.parse(varargin{:});

      % Get or estimate Nmax from the inputs
      if isempty(p.Results.Nmax)
        k_medium = ott.Tmatrix.parser_k_medium(p, 2*pi);
        Nmax = ott.utils.ka2nmax(shape.maxRadius * k_medium);
      else
        Nmax = p.Results.Nmax;

        % We only support square matricies for now
        if numel(Nmax) ~= 1
          Nmax = max(Nmax(:));
        end
      end

      % Get the symmetry of the shape
      [~,~, z_rotational_symmetry] = shape.axialSymmetry();
      [~,~, z_mirror_symmetry] = shape.mirrorSymmetry();

      % Get the coordinates of the shape
      rtp = shape.angulargrid(Nmax);
      normals = shape.normals(rtp(:, 2), rtp(:, 3));

      % inputParser will take the last parameter, so varargin just needs to
      % be before the replacements for varargin.
      tmatrix = ott.TmatrixPm(rtp, normals, varargin{:}, ...
          'Nmax', Nmax, ...
          'z_mirror_symmetry', z_mirror_symmetry, ...
          'z_rotational_symmetry', z_rotational_symmetry);
    end
  end

  methods (Access=protected)

    function [coeff_matrix, incident_wave_matrix] = setup(tmatrix, ...
        Nmax, rtp, normals, progress)
      % SETUP calculates the coefficient and incident wave matrices

      npoints = size(rtp, 1);
      total_orders = ott.utils.combined_index(Nmax, Nmax);

      r = rtp(:, 1);
      theta = rtp(:, 2);
      phi = rtp(:, 3);

      % 3 vector components at each point, c/d,p/q coefficient per order
      coeff_matrix = zeros(6*npoints,4*total_orders);
      incident_wave_matrix = zeros(6*npoints,2*total_orders);

      k_relative = tmatrix.k_particle/tmatrix.k_medium;

      import ott.utils.vswf;
      import ott.utils.perpcomponent;

      for n = 1:Nmax
        for m = -n:n

          % INCIDENT-SCATTERED
          [M1,N1,~,~,M2,N2] = vswf(n,m,tmatrix.k_medium*r,theta,phi);
          [M3,N3] = vswf(n,m,tmatrix.k_particle*r,theta,phi,3);

          ci = ott.utils.combined_index(n,m);

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

          % Output progress
          if progress ~= 0 && mod(ci, progress) == 0
            disp(['TmatrixPM:setup (' num2str(ci) '/' ...
                num2str(total_orders) ') ' ...
                num2str(floor(ci/total_orders*100.0)) '%']);
          end
        end
      end

    end

  end

  methods
    function tmatrix = TmatrixPm(rtp, normals, varargin)
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
      %  TMATRIXPM(..., 'z_rotational_symmetry', sym) if true the particle
      %  is assumed to be rotationally symmetric about the z-axis and
      %  phi must be all the same angle.
      %
      %  TMATRIXPM(..., 'z_mirror_symmetry', sym) not yet implemented.

      tmatrix = tmatrix@ott.Tmatrix();

      % TODO: Different T and T2 size
      % TODO: Different row/column Nmax

      % Parse inputs
      p = inputParser;
      p.addParameter('Nmax', []);
      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('k_particle', []);
      p.addParameter('wavelength_particle', []);
      p.addParameter('index_particle', []);
      p.addParameter('index_relative', []);
      p.addParameter('wavelength0', []);
      p.addParameter('internal', false);
      p.addParameter('z_rotational_symmetry', false);
      p.addParameter('z_mirror_symmetry', false);
      p.addParameter('progress', []);
      p.parse(varargin{:});

      % Store inputs k_medium and k_particle
      [tmatrix.k_medium, tmatrix.k_particle] = ...
          tmatrix.parser_wavenumber(p, 2*pi);

      % Get or estimate Nmax from the inputs
      if isempty(p.Results.Nmax)
        Nmax = ka2nmax(max(rtp(:, 1)) * tmatrix.k_medium);
      else
        Nmax = p.Results.Nmax;

        % We only support square matricies for now
        if numel(Nmax) ~= 1
          Nmax = max(Nmax(:));
        end
      end

      % Parse progress input
      if ~isempty(p.Results.progress)
        progress = p.Results.progress;
        progress(progress < 1) = floor(...
            progress(progress < 1) * ott.utils.combined_index(Nmax, Nmax));
        if numel(progress) == 1
          progress = [progress, progress];
        end
      else
        progress = [0, 0];
      end

      % Calculate coefficient and incident wave matrices
      [coeff_matrix, incident_wave_matrix] = tmatrix.setup(...
          Nmax, rtp, normals, progress(1));

      total_orders = ott.utils.combined_index(Nmax, Nmax);

      import ott.utils.*

      % Generate T-matrix

      if p.Results.z_rotational_symmetry == 0

        % Infinite rotational symmetry

        T = sparse(2*total_orders,2*total_orders);
        T2 = sparse(2*total_orders,2*total_orders);

        [n, m] = ott.utils.combined_index((1:total_orders).');

        for mi = -Nmax:Nmax

          % Scattered modes have the same m as incident modes
          axial_modes = m == mi;
          modes = [ axial_modes; axial_modes ];

          % This is for future, using different T and T2 size
          emodes = modes;
          omodes = modes;
          iomodes = [ omodes; emodes ];

          if p.Results.z_mirror_symmetry

            % Correct the incident modes to include even/odd modes
            even_modes = logical(mod(n + mi, 2));
            imodes_evn = modes & [ even_modes; ~even_modes ];
            imodes_odd = modes & [ ~even_modes; even_modes ];

            % Correct the outgoing modes to include even/odd modes
            even_modes = logical(mod(n + m, 2));
            omodes_evn = omodes & [ even_modes; ~even_modes ];
            omodes_odd = omodes & [ ~even_modes; even_modes ];
            emodes_evn = emodes & [ even_modes; ~even_modes ];
            emodes_odd = emodes & [ ~even_modes; even_modes ];
            iomodes_evn = [ omodes_evn; emodes_evn ];
            iomodes_odd = [ omodes_odd; emodes_odd ];

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_evn);
            Tcol = coeff_matrix(:, iomodes_evn) \ incident_wave_vectors;
            T(omodes_evn, imodes_evn) = Tcol(1:sum(omodes_evn), :);
            T2(emodes_evn, imodes_evn) = Tcol((1+sum(omodes_evn)):end, :);

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_odd);
            Tcol = coeff_matrix(:, iomodes_odd) \ incident_wave_vectors;
            T(omodes_odd, imodes_odd) = Tcol(1:sum(omodes_odd), :);
            T2(emodes_odd, imodes_odd) = Tcol((1+sum(omodes_odd)):end, :);

          else

            % Scatter the modes
            incident_wave_vectors = incident_wave_matrix(:, modes);
            Tcol = coeff_matrix(:, iomodes) \ incident_wave_vectors;
            T(omodes, modes) = Tcol(1:sum(omodes), :);
            T2(emodes, modes) = Tcol((1+sum(omodes)):end, :);

          end

        end

      elseif p.Results.z_rotational_symmetry ~= 1

        % Discrete rotational symmetry

        T = sparse(2*total_orders,2*total_orders);
        T2 = sparse(2*total_orders,2*total_orders);

        [n, m] = ott.utils.combined_index((1:total_orders).');

        for mi = -Nmax:Nmax

          % Calculate which modes preseve symmetry, m = +/- ip
          axial_modes = mod(m - mi, p.Results.z_rotational_symmetry) == 0;
          incm_modes = m == mi;
          modes = [ axial_modes; axial_modes ];
          imodes = [ incm_modes; incm_modes ];

          % This is for future, using different T and T2 size
          emodes = modes;
          omodes = modes;
          iomodes = [ omodes; emodes ];

          if p.Results.z_mirror_symmetry

            % Correct the incident modes to include even/odd modes
            even_modes = logical(mod(n + mi, 2));
            imodes_evn = modes & [ even_modes; ~even_modes ];
            imodes_odd = modes & [ ~even_modes; even_modes ];

            % Correct the outgoing modes to include even/odd modes
            even_modes = logical(mod(n + m, 2));
            omodes_evn = omodes & [ even_modes; ~even_modes ];
            omodes_odd = omodes & [ ~even_modes; even_modes ];
            emodes_evn = emodes & [ even_modes; ~even_modes ];
            emodes_odd = emodes & [ ~even_modes; even_modes ];
            iomodes_evn = [ omodes_evn; emodes_evn ];
            iomodes_odd = [ omodes_odd; emodes_odd ];

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_evn);
            Tcol = coeff_matrix(:, iomodes_evn) \ incident_wave_vectors;
            T(omodes_evn, imodes_evn) = Tcol(1:sum(omodes_evn), :);
            T2(emodes_evn, imodes_evn) = Tcol((1+sum(omodes_evn)):end, :);

            % Solve for the even scattered modes
            incident_wave_vectors = incident_wave_matrix(:, imodes_odd);
            Tcol = coeff_matrix(:, iomodes_odd) \ incident_wave_vectors;
            T(omodes_odd, imodes_odd) = Tcol(1:sum(omodes_odd), :);
            T2(emodes_odd, imodes_odd) = Tcol((1+sum(omodes_odd)):end, :);

          else

            % Scatter the modes
            incident_wave_vectors = incident_wave_matrix(:, imodes);
            Tcol = coeff_matrix(:, iomodes) \ incident_wave_vectors;
            T(omodes, imodes) = Tcol(1:sum(omodes), :);
            T2(emodes, imodes) = Tcol((1+sum(omodes)):end, :);

          end

        end

      else

        if p.Results.z_mirror_symmetry

          % Only mirror symmetry
          % Parity is conserved, even modes go to even modes, etc.
          % Reference: https://arxiv.org/pdf/physics/0702045.pdf

          T = sparse(2*total_orders,2*total_orders);
          T2 = sparse(2*total_orders,2*total_orders);

          [n, m] = ott.utils.combined_index((1:total_orders).');
          even_modes = logical(mod(n + m, 2));
          modes = [ even_modes; ~even_modes ];

          % This is for future, using different T and T2 size
          imodes = modes;
          omodes = modes;
          iomodes = [ omodes; imodes ];

          % Solve for the even scattered modes
          incident_wave_vectors = incident_wave_matrix(:, modes);
          Tcol = coeff_matrix(:, iomodes) \ incident_wave_vectors;
          T(omodes, modes) = Tcol(1:sum(omodes), :);
          T2(imodes, modes) = Tcol((1+sum(omodes)):end, :);

          % Solve for the odd scattered modes
          incident_wave_vectors = incident_wave_matrix(:, ~modes);
          Tcol = coeff_matrix(:, ~iomodes) \ incident_wave_vectors;
          T(~omodes, ~modes) = Tcol(1:sum(~omodes),:);
          T2(~imodes, ~modes) = Tcol((1+sum(~omodes)):end,:);

        else

          % No rotational or mirror symmetry

          Tcol = coeff_matrix \ incident_wave_matrix;
          T = Tcol(1:2*total_orders,:);
          T2 = Tcol((1+2*total_orders):4*total_orders,:);

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

      % Store the type of T-matrix
      tmatrix.type = 'scattered';
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
