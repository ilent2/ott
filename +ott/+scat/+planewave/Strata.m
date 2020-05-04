classdef Strata
% Describes scattering of a plane wave by layered (stratified) planes
%
% Implementation of the method described in
%
%   James R. Wait, Transmission and Reflection of Electromagnetic
%   Waves in the Presence of Stratified Media.
%   Journal of Research of the National Bureau of Standards
%   Vol. 61, No.3, September 1958 Research Paper 2899
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    normal          % Vector perpendicular to the surfaces
    index_relative  % Relative refractive index of each layer
  end

  methods (Static)
    function [aS, bS] = coeffS(normal, depths, kvecs)
      % Calculate the S-scattering coefficients (out-of-plane)

      % For E-H relationship (and units!!!)
      sig = 1.0;

      num_ir = numel(depths) + 1;
      
      % Calcuate term related to wave-vector
%       unormal = sqrt(-dot(kvecs, kvecs) + vecnorm(cross(kvecs(:, 1), normal)).^2);
      knormal = dot(kvecs, repmat(normal, 1, num_ir + 1));
      unormal = 1i.*knormal;

      % Calculate terms for boundary matrix
      M = exp(-unormal .* [0, 0, depths]);
      N = exp(unormal .* [0, 0, depths]);

      % Calculate transverse field terms (for boundary matrix)
      Mt = -1i ./ sig .* unormal .* M;
      Nt = 1i ./ sig .* unormal .* N;

      % Calculate transfer terms between layers
      dd = [1, exp(-1i.*depths.*knormal(:, 2:end-1))];

      % Assemble a matrix of boundary conditions between the layers
      A = zeros(2*num_ir, num_ir+3);
      for ii = 1:num_ir
        col_idx = 2*(ii-1)+1;
        A(ii, col_idx:col_idx+3) = ...
            [M(ii).*dd(ii), N(ii).*dd(ii), -M(ii+1), -N(ii+1)];
        A(ii+num_ir, col_idx:col_idx+3) = ...
            [Mt(ii).*dd(ii), Nt(ii).*dd(ii), -Mt(ii+1), -Nt(ii+1)];
      end

      % Reshape A so its in the form A x = b
      % First column of A is our ab vector (except last and first element)
      % Remaining columns of A are A matrix

      b = -A(:, 1);
      A = A(:, 2:end-1);      % Discard: a0 = 1, bM = 0

      % Solve linear over-determined system
      ab = A \ b;

      % Pack results into aS and bS
      aS = ab(2:2:end);
      bS = ab(1:2:end);
    end

    function [aP, bP] = coeffP(normal, depths, kvecs)
      % Calculate the P-scattering coefficients (in-plane)

      % For E-H relationship (and units!!!)
      n_rel = vecnorm(kvecs) ./ vecnorm(kvecs(:, 1));
      sig = n_rel.^2;
%       sig = 1.0;

      num_ir = numel(depths) + 1;
      
      % Calcuate term related to wave-vector
%       unormal = sqrt(-dot(kvecs, kvecs) + vecnorm(cross(kvecs(:, 1), normal)).^2);
      knormal = dot(kvecs, repmat(normal, 1, num_ir + 1));
      unormal = 1i.*knormal;

      % Calculate terms for boundary matrix
      % Based on Wait Eq. 2
      M = exp(-unormal .* [0, 0, depths]);
      N = exp(unormal .* [0, 0, depths]);

      % Calculate transverse field terms (for boundary matrix)
      Mt = -1i ./ sig .* unormal .* M;
      Nt = 1i ./ sig .* unormal .* N;

      % Calculate transfer terms between layers
      dd = [1, exp(-1i.*depths.*knormal(:, 2:end-1))];

      % Assemble a matrix of boundary conditions between the layers
      % Based on Wait Eq. 4
      A = zeros(2*num_ir, num_ir+3);
      for ii = 1:num_ir
        col_idx = 2*(ii-1)+1;
        A(ii, col_idx:col_idx+3) = ...
            [M(ii).*dd(ii), N(ii).*dd(ii), -M(ii+1), -N(ii+1)];
        A(ii+num_ir, col_idx:col_idx+3) = ...
            [Mt(ii).*dd(ii), Nt(ii).*dd(ii), -Mt(ii+1), -Nt(ii+1)];
      end

      % Reshape A so its in the form A x = b
      % First column of A is our ab vector (except last and first element)
      % Remaining columns of A are A matrix

      b = -A(:, 1);
      A = A(:, 2:end-1);      % Discard: a0 = 1, bM = 0

      % Solve linear over-determined system
      ab = A \ b;

      % Pack results into aS and bS (and convert from H to E)
      aP = ab(2:2:end) ./ n_rel(2:end).';
      bP = ab(1:2:end) ./ n_rel(1:end-1).';
    end
  end

  methods
    function strata = Strata(normal, index_relative)
      strata.normal = normal;
      strata.index_relative = index_relative;
    end

    function [rbeam, tbeam, ibeams] = scatter(strata, beam)
      % Calculate the reflected, transmitted and internal beams
      %
      % Usage
      %   [rbeam, tbeam, ibeams] = strata.scatter(beam)

      % Cast the beam to a plane wave
      if ~isa(beam, 'ott.beam.abstract.PlaneWave')
        beam = ott.beam.abstract.PlaneWave(beam);
      end

      % Get incident wave-vector (direction)
      ki = beam.wavevector;

      % Calculate reflected wave-vector (direction)
      % wavenumber is the same since the medium is the same
      % Direction of the normal component changes (reflected)
      kr = ki - 2.*dot(plane.normal, ki).*plane.normal;

      % Calculate wave vectors in each layer
      num_ir = numel(strata.index_relative);
      kt = zeros(3, num_ir);
      ir_prev = 1.0;
      k_prev = ki;
      for ii = 1:numel(num_ir);

        % Calculate relative index between two mediums
        index_relative = strata.index_relative(ii)./ir_prev;
        ir_prev = strata.index_relative(ii);

        % Orthogonal components are unchanged
        % Normal component is scaled by relative index
        ko = k_prev - dot(strata.normal, k_prev).*strata.normal; % orthogonal
        kn2 = index_relative.^2 .* dot(k_prev, k_prev) - dot(ko, ko); % normal
        kt(:, ii) = ko + sqrt(kn2) .* strata.normal;
      end

      % Solve the boundary value problem
      [aS, bS] = strata.coeffS(strata.normal, strata.depth, [ki, kt]);
      [aP, bP] = strata.coeffP(strata.normal, strata.depth, [ki, kt]);

      polarisation = beam.polarisation;

      % TODO: Use field vector

      % Split the field (polarisation) into s and p vectors
      svec = cross(beam.direction, plane.normal);
      Es = dot(polarisation, svec);
      Ep = vecnorm(polarisation - Es .* svec);

      % Calculate transmitted and reflected pvec directions
      pvecr = cross(svec, kr)./vecnorm(kr);
      pvect = cross(svec, kt)./vecnorm(kt);

      rbeam = ott.beam.abstract.PlaneWave('wave_vector', kr, ...
          'polarisation', bS(1) .* Es .* svec + bP(1) .* Ep .* pvecr, ...
          'index_medium', beam.index_medium, ...
          'origin', plane.position);
      tbeam = ott.beam.abstract.PlaneWave('wave_vector', kt, ...
          'polarisation', aS(end) .* Es .* svec + aP(end) .* Ep .* pvect, ...
          'index_medium', beam.index_medium .* plane.index_relative, ...
          'origin', plane.position);

      if nargout == 3

        ibeams(2*(length(aS)-1)) = ott.beam.abstract.PlaneWave();
        for ii = 2:length(aS)
          ibeams(ii-1) = ott.beam.abstract.PlaneWave(...
            'wave_vector', kt(:, ii-1), ...
            'polarisation', aS(ii-1) .* Es .* svec + aP(ii-1) .* Ep .* pvect, ...
            'index_medium', beam.index_medium .* strata.index_relative(ii), ...
            'origin', beam.origin);

          ibeams(ii+length(aS)-2) = ott.beam.abstract.PlaneWave(...
            'wave_vector', kt(:, ii-1), ...
            'polarisation', bS(ii) .* Es .* svec + bP(ii) .* Ep .* pvect, ...
            'index_medium', beam.index_medium .* strata.index_relative(ii), ...
            'origin', beam.origin);
        end
      end
    end
  end
end
