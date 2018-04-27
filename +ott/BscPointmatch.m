classdef BscPointmatch < ott.Bsc
%BscPointmatch base class for BSC generated using point matching
% Provides support for both farfield and focal plane point matching.
%
% BscPointmatch properties:
%
% BscPointmatch static methods:
%   bsc_farfield          Does point matching in the farfield
%   bsc_focalplane        Does point matching around the focal plane
%
% Based on bsc_pointmatch_focalplane and bsc_pointmatch_farfield
% from version 1 of the optical tweezers toolbox.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Static, Access=protected)

    function [a, b] = bsc_farfield(nn, mm, e_field, theta, phi, ...
        zero_rejection_level)
      % FARFIELD point match beam coefficients in farfield
      %
      % nn, mm are the mode indices to include in the coefficient matrix.
      %
      % e_field is the E-field to point match.
      %
      % theta, phi are the coordinates of the Efield values.
      %
      % zero_rejection_level (optional) removes modes with less
      % power than zero_rejection_level.  Default 1e-8.

      if nargin < 6
        zero_rejection_level = 1e-8;
      end

      % Generate coefficient matrix
      coefficient_matrix = zeros(length(e_field), 2*length(nn));
      for n = 1:max(nn)
        ci=find(nn==n);

        [~,dtY,dpY]= ott.utils.spharm(n,mm(ci),theta,phi);

        coefficient_matrix(:,ci) = [dpY;-dtY] * 1i^(n+1)/sqrt(n*(n+1));
        coefficient_matrix(:,ci+length(nn)) = [dtY;dpY] * 1i^(n)/sqrt(n*(n+1));
      end

      % Do point matching
      expansion_coefficients = coefficient_matrix \ e_field;

      % Unpack results into a and b vectors
      fa = expansion_coefficients(1:end/2,:);
      fb = expansion_coefficients(1+end/2:end,:);

      % Look for non-zero elements, only keep non-zeros
      pwr = abs(fa).^2+abs(fb).^2;
      non_zero = pwr>zero_rejection_level*max(pwr);
      nn=nn(non_zero);
      mm=mm(non_zero);
      fa=fa(non_zero);
      fb=fb(non_zero);

      % Make the beam vector and store the coefficients
      [a, b] = ott.Bsc.make_beam_vector(fa, fb, nn, mm);
    end

    function [a, b] = bsc_focalplane(nn, mm, e_field, kr, theta, phi, ...
        zero_rejection_level)
      % FOCALPLANE point match beam coefficients around focal plane
      %
      % nn, mm are the mode indices to include in the coefficient matrix.
      %
      % e_field is the E-field to point match.
      %
      % kr, theta, phi are the coordinates of the Efield values.
      %
      % zero_rejection_level (optional) removes modes with less
      % power than zero_rejection_level.  Default 1e-8.

      if nargin < 7
        zero_rejection_level = 1e-8;
      end

      % Generate coefficient matrix
      coefficient_matrix = zeros(length(e_field), length(nn));
      for n = 1:length(nn)

         % Find RgM, RgN as appropriate for each mode
         [M,N] = vswfcart(nn(n),mm(n),kr,theta,phi,3);
         if rem(nn(n),2) == 0
            % Even n
            MN = [ M(:,1); M(:,2); M(:,3) ];
         else
            % Odd n
            MN = [ N(:,1); N(:,2); N(:,3) ];
         end
         coefficient_matrix(:,n) = MN;

      end

      % Solve the linear system
      expansion_coefficients = coefficient_matrix \ e_field;

      % Look for non-zero elements, only keep non-zeros
      non_zero = abs(expansion_coefficients) ...
          > max(abs(expansion_coefficients)) * zero_rejection_level;
      expansion_coefficients = expansion_coefficients(non_zero);
      nn = nn(non_zero);
      mm = mm(non_zero);

      % Calculate beam vectors
      fa = zeros(size(nn));
      fb = zeros(size(nn));
      for n = 1:length(nn)

         if rem(nn(n),2) == 0
            fa(n) = expansion_coefficients(n);
            fb(n) = expansion_coefficients(n) * sign(mm(n));
         else
            fa(n) = expansion_coefficients(n) * sign(mm(n));
            fb(n) = expansion_coefficients(n);
         end

      end

      % Make the beam vector and store the coefficients
      [a, b] = ott.Bsc.make_beam_vector(fa, fb, nn, mm);
    end
  end

  methods
    function beam = BscPointmatch(varargin)
      beam = beam@ott.Bsc();
    end
  end
end
