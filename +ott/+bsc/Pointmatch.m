classdef BscPointmatch < ott.Bsc
%BscPointmatch base class for BSC generated using point matching
% Provides support for both far-field and focal plane point matching.
%
% Properties
%   inv_coefficient_matrix  Pseudo-inverse coefficient matrix for PM
%   a            (Bsc) Beam shape coefficients a vector
%   b            (Bsc) Beam shape coefficients b vector
%   type         (Bsc) Beam type (incident, scattered, total)
%   basis        (Bsc) VSWF beam basis (incoming, outgoing or regular)
%   Nmax         (Bsc) Truncation number for VSWF coefficients
%   power        (Bsc) Power of the beam [M*L^2/S^3]
%   Nbeams       (Bsc) Number of beams in this Bsc object
%   wavelength   (Bsc) Wavelength of beam [L]
%   speed        (Bsc) Speed of beam in medium [L/T]
%   omega        (Bsc) Angular frequency of beam [2*pi/T]
%   k_medium     (Bsc) Wavenumber in medium [2*pi/L]
%   dz           (Bsc) Absolute cumulative distance the beam has moved
%
% Methods
%   cleanCoefficientMatrix  Removes coefficient matrix data from the beam.
%
% Static methods
%   bsc_farfield          Does point matching in the farfield
%   bsc_focalplane        Does point matching around the focal plane
%
% See also bsc_farfield, bsc_focalplane and ott.BscPmGauss.
%
% Based on bsc_pointmatch_focalplane and bsc_pointmatch_farfield
% from version 1 of the optical tweezers toolbox.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Static)

    function [a, b, icm] = bsc_farfield(nn, mm, e_field, theta, phi, ...
        varargin)
      % point match beam coefficients in farfield
      %
      % [a, b, icm] = bsc_farfield(nn, mm, e_field, theta, phi, ...)
      % a, b are the beam coefficients.  cm is the coefficient matrix.
      %
      % nn, mm are the mode indices to include in the coefficient matrix.
      %
      % e_field is the E-field to point match.
      % The format should be [ Etheta(:); Ephi(:) ]
      %
      % theta, phi are the coordinates of the Efield values.
      %
      % Optional named arguments:
      %   zero_rejection_level   val   removes modes with less power than
      %       zero_rejection_level.  Default 1e-8.  Use [] to disable.
      %   inv_coefficient_matrix     mat   Coefficient matrix to use.
      %       default [].
      %   invert_coefficient_matrix  bool  True to invert coefficient
      %       matrix for point matching.  Default nargout == 3.

      p = inputParser;
      p.addParameter('inv_coefficient_matrix', []);
      p.addParameter('zero_rejection_level', 1e-8);
      p.addParameter('invert_coefficient_matrix', nargout == 3);
      p.parse(varargin{:});

      % Generate coefficient matrix
      icm = p.Results.inv_coefficient_matrix;
      if isempty(icm)
        coefficient_matrix = zeros(length(e_field), 2*length(nn));
        for n = 1:max(nn)
          ci=find(nn==n);

          [~,dtY,dpY]= ott.utils.spharm(n,mm(ci),theta,phi);

          coefficient_matrix(:,ci) = [dpY;-dtY] * 1i^(n+1)/sqrt(n*(n+1));
          coefficient_matrix(:,ci+length(nn)) = [dtY;dpY]*1i^(n)/sqrt(n*(n+1));
        end

        % Invert coefficient matrix for icm
        if nargout == 3 || p.Results.invert_coefficient_matrix
          icm = pinv(coefficient_matrix);
        end

        % Do point matching
        if p.Results.invert_coefficient_matrix
          expansion_coefficients = icm * e_field;
        else
          expansion_coefficients = coefficient_matrix \ e_field;
        end

      else
        assert(size(icm, 2) == length(e_field), ...
            'Number of cols in coefficient matrix must match length(e_field)');

        % Do point matching
        expansion_coefficients = icm * e_field;
      end

      % Unpack results into a and b vectors
      fa = expansion_coefficients(1:end/2,:);
      fb = expansion_coefficients(1+end/2:end,:);

      % Look for non-zero elements, only keep non-zeros
      if ~isempty(p.Results.zero_rejection_level)
        pwr = abs(fa).^2+abs(fb).^2;
        non_zero = pwr>p.Results.zero_rejection_level*max(pwr);
        nn=nn(non_zero);
        mm=mm(non_zero);
        fa=fa(non_zero);
        fb=fb(non_zero);
      end

      % Make the beam vector and store the coefficients
      [a, b] = ott.Bsc.make_beam_vector(fa, fb, nn, mm);
    end

    function [a, b, cm] = bsc_focalplane(nn, mm, e_field, kr, theta, phi, ...
        varargin)
      % point match beam coefficients around focal plane
      %
      % [a, b, cm] = bsc_focalplane(nn, mm, e_field, kr, theta, phi, ...)
      % a, b are the beam coefficients.  cm is the coefficient matrix.
      %
      % nn, mm are the mode indices to include in the coefficient matrix.
      %
      % e_field is a vector of E-field to points to match.
      % The format should be [ Ex(:); Ey(:); Ez(:) ]
      %
      % kr, theta, phi are the coordinates of the Efield values.
      % These should be vectors of the same as length(e_field)/3.
      %
      % Optional named arguments:
      %   zero_rejection_level   val   removes modes with less power than
      %       zero_rejection_level.  Default 1e-8.  Use [] to disable.
      %   coefficient_matrix     mat   Coefficient matrix to use.
      %       default [].

      p = inputParser;
      p.addParameter('inv_coefficient_matrix', []);
      p.addParameter('zero_rejection_level', 1e-8);
      p.addParameter('invert_coefficient_matrix', nargout == 3);
      p.parse(varargin{:});
      
      assert(length(e_field) == numel(e_field), ...
        'e_field must be N element vector');
      assert(numel(e_field)/3 == numel(kr), 'kr must be same size as e_field/3');
      assert(numel(e_field)/3 == numel(theta), 'theta must be same size as e_field/3');
      assert(numel(e_field)/3 == numel(phi), 'phi must be same size as e_field/3');

      % Generate coefficient matrix
      icm = p.Results.inv_coefficient_matrix;
      if isempty(icm)
        coefficient_matrix = zeros(length(e_field), length(nn));
        for n = 1:length(nn)

           % Find RgM, RgN as appropriate for each mode
           [M,N] = ott.utils.vswfcart(nn(n),mm(n),kr,theta,phi,3);
           if rem(nn(n),2) == 0
              % Even n
              MN = [ M(:,1); M(:,2); M(:,3) ];
           else
              % Odd n
              MN = [ N(:,1); N(:,2); N(:,3) ];
           end
           coefficient_matrix(:,n) = MN;

        end

        % Invert coefficient matrix for icm
        if nargout == 3 || p.Results.invert_coefficient_matrix
          icm = pinv(coefficient_matrix);
        end

        % Do point matching
        if p.Results.invert_coefficient_matrix
          expansion_coefficients = icm * e_field;
        else
          expansion_coefficients = coefficient_matrix \ e_field;
        end
      else
        assert(size(icm, 2) == length(e_field), ...
            'Number of rows in coefficient matrix must match length(e_field)');

        % Do point matching
        expansion_coefficients = icm * e_field;
      end

      % Look for non-zero elements, only keep non-zeros
      if ~isempty(p.Results.zero_rejection_level)
        non_zero = abs(expansion_coefficients) ...
            > max(abs(expansion_coefficients)) * p.Results.zero_rejection_level;
        expansion_coefficients = expansion_coefficients(non_zero);
        nn = nn(non_zero);
        mm = mm(non_zero);
      end

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

  properties (SetAccess=protected)
    inv_coefficient_matrix      % Coefficient matrix used in point matching
  end

  methods (Access=protected)
    function beam = BscPointmatch(varargin)
      % Protected constructor for BscPointmatch object
      %
      % See also ott.BscPmGauss and ott.BscPmParaxial
      beam = beam@ott.Bsc();
    end
  end

  methods
    function cleanCoefficientMatrix(beam)
      % Remove the coefficient matrix data
      beam.inv_coefficient_matrix = [];
    end
  end
end
