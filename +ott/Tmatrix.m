classdef Tmatrix < matlab.mixin.Copyable
%Tmatrix abstract class representing T-matrix of a scattering particle
%
% Tmatrix properties:
%   data          The T-matrix this class encapculates
%
% Tmatrix methods:
%
% Static methods:
%
% Abstract methods:
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

 properties (SetAccess=protected)
  data            % The matrix this class encapsulates
 end

  properties (Dependent)
    Nmax          % Current size of T-matrix
  end

  methods (Static)
    function k_medium = parser_k_medium(p)
      %PARSER_K_MEDIUM helper to get k_medium from a parser object

      if ~isempty(p.Results.k_medium)
        k_medium = p.Results.k_medium;
      elseif ~isempty(p.Results.wavelength_medium)
        k_medium = 2.0*pi/p.Results.wavelength_medium;
      elseif ~isempty(p.Results.index_medium)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_medium');
        end
        k_medium = p.Results.index_medium*2.0*pi/p.Results.wavelength0;
      else
        error('Unable to determine k_medium from inputs');
      end
    end

    function k_particle = parser_k_particle(p)
      %PARSER_K_PARTICLE helper to get k_particle from a parser object

      if ~isempty(p.Results.k_particle)
        k_particle = p.Results.k_particle;
      elseif ~isempty(p.Results.wavelength_particle)
        k_particle = 2.0*pi/p.Results.wavelength_particle;
      elseif ~isempty(p.Results.index_particle)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_particle');
        end
        k_particle = p.Results.index_particle ...
            * 2.0*pi/p.Results.wavelength0;
      else
        error('Unable to determine k_particle from inputs');
      end
    end
  end

 methods (Abstract)
 end

 methods (Access=protected)
  function tmatrix = Tmatrix(data)
    %TMATRIX construct a new T-matrix object
    %
    % TMATRIX() leaves the data uninitialised.
    %
    % TMATRIX(data) initializes the data with the matrix data.

    if nargin == 1
      tmatrix.data = data;
    end
  end
 end

  methods
    function nmax = get.Nmax(tmatrix)
      %get.Nmax calculate Nmax from the current T-matrix data
      nmax1 = ott.utils.combined_index(size(tmatrix.data, 1));
      nmax2 = ott.utils.combined_index(size(tmatrix.data, 2));

      % Support non-square T-matrices
      if nmax1 == nmax2
        nmax = nmax1;
      else
        nmax = [nmax1 nmax2];
      end
    end

    function tmatrix = set.Nmax(tmatrix, nmax)
      %set.Nmax resizes the T-matrix

      if length(nmax) == 2
        nmax1 = nmax(1);
        nmax2 = nmax(2);
      else
        nmax1 = nmax(1);
        nmax2 = nmax(1);
      end

      total_orders1 = ott.utils.combined_index(nmax1, nmax1);
      total_orders2 = ott.utils.combined_index(nmax2, nmax2);

      midpoint1 = size(tmatrix.data, 1)/2;
      midpoint2 = size(tmatrix.data, 2)/2;

      % Split T-matrix into quadrants
      A11 = tmatrix.data(1:midpoint1, 1:midpoint2);
      A12 = tmatrix.data(1:midpoint1, (midpoint2+1):end);
      A21 = tmatrix.data((midpoint1+1):end, 1:midpoint2);
      A22 = tmatrix.data((midpoint1+1):end, (midpoint2+1):end);

      % Resize rows
      if total_orders1 > midpoint1

        [row_index,col_index,a] = find(A11);
        A11 = sparse(row_index,col_index,a,total_orders1,midpoint2);

        [row_index,col_index,a] = find(A12);
        A12 = sparse(row_index,col_index,a,total_orders1,midpoint2);

        [row_index,col_index,a] = find(A21);
        A21 = sparse(row_index,col_index,a,total_orders1,midpoint2);

        [row_index,col_index,a] = find(A22);
        A22 = sparse(row_index,col_index,a,total_orders1,midpoint2);

      else

        A11 = A11(1:total_orders1, :);
        A12 = A12(1:total_orders1, :);
        A21 = A21(1:total_orders1, :);
        A22 = A22(1:total_orders1, :);

      end

      % Resize cols
      if total_orders2 > midpoint2

        [row_index,col_index,a] = find(A11);
        A11 = sparse(row_index,col_index,a,total_orders1,total_orders2);

        [row_index,col_index,a] = find(A12);
        A12 = sparse(row_index,col_index,a,total_orders1,total_orders2);

        [row_index,col_index,a] = find(A21);
        A21 = sparse(row_index,col_index,a,total_orders1,total_orders2);

        [row_index,col_index,a] = find(A22);
        A22 = sparse(row_index,col_index,a,total_orders1,total_orders2);

      else

        A11 = A11(:, 1:total_orders2);
        A12 = A12(:, 1:total_orders2);
        A21 = A21(:, 1:total_orders2);
        A22 = A22(:, 1:total_orders2);

      end

      magA = full(sum(sum(abs(tmatrix.data).^2)));

      % Recombined T-matrix from quadrants
      tmatrix.data = [ A11 A12; A21 A22 ];

      magB = full(sum(sum(abs(tmatrix.data).^2)));
      apparent_error = abs( magA - magB )/magA;

      warning_error_level = 1e-6;
      if apparent_error > warning_error_level
          warning('ott:Tmatrix:set.Nmax:truncation', ...
              ['Apparent error of ' num2str(apparent_error)]);
      end
    end
  end
end
