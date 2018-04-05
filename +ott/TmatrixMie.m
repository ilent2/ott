classdef TmatrixMie < ott.Tmatrix
%TmatrixMie construct T-matrix from Mie scattering coefficients
%
% TmatrixMie properties:
%   radius            The radius of the sphere the T-matrix represents
%   k_medium          Wavenumber in the trapping medium
%   k_particle        Wavenumber of the particle
%
% This class is based on tmatrix_mie.m from ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    radius            % Radius of particle
    k_medium          % Wavenumber of medium
    k_particle        % Wavenumber of particle
  end

  methods
    function tmatrix = TmatrixMie(radius, varargin)
      %TMATRIXMIE construct a new Mie T-matrix for a sphere with size radius.
      %
      %  TMATRIXMIE(..., 'Nmax', Nmax) specifies the size of the
      %  T-matrix to use.  If not specified, the size is calculated
      %  from ott.utils.ka2nmax(radius*k_medium).
      %
      %  TMATRIXMIE(..., 'k_medium', k)
      %  or TMATRIXMIE(..., 'wavelength_medium', wavelength)
      %  or TMATRIXMIE(..., 'index_medium', index)
      %  specify the wavenumber, wavelength or index in the medium.
      %
      %  TMATRIXMIE(..., 'k_particle', k)
      %  or TMATRIXMIE(..., 'wavelength_particle', wavelength)
      %  or TMATRIXMIE(..., 'index_particle', index)
      %  specify the wavenumber, wavelength or index in the particle.
      %
      %  TMATRIXMIE(..., 'wavelength0', wavelength) specifies the
      %  wavelength in the vecuum, required when index_particle or
      %  index_medium are specified.
      %
      %  TMATRIXMIE(..., 'internal', internal) if true, calculates the
      %  T-matrix for the internal coefficients.

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
      p.parse(varargin{:});

      % Store inputs: radius
      tmatrix.radius = radius;

      % Store inputs: k_medium
      if ~isempty(p.Results.k_medium)
        tmatrix.k_medium = p.Results.k_medium;
      elseif ~isempty(p.Results.wavelength_medium)
        tmatrix.k_medium = 2.0*pi/p.Results.wavelength_medium;
      elseif ~isempty(p.Results.index_medium)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_medium');
        end
        tmatrix.k_medium = p.Results.index_medium*2.0*pi/p.Results.wavelength0;
      else
        error('Unable to determine k_medium from inputs');
      end

      % Store inputs: k_particle
      if ~isempty(p.Results.k_particle)
        tmatrix.k_particle = p.Results.k_particle;
      elseif ~isempty(p.Results.wavelength_particle)
        tmatrix.k_particle = 2.0*pi/p.Results.wavelength_particle;
      elseif ~isempty(p.Results.index_particle)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_particle');
        end
        tmatrix.k_particle = p.Results.index_particle ...
            * 2.0*pi/p.Results.wavelength0;
      else
        error('Unable to determine k_particle from inputs');
      end

      % If Nmax not specified, choose a good Nmax
      if isempty(p.Results.Nmax)
        Nmax = ott.utils.ka2nmax(tmatrix.k_medium*radius);
      else
        Nmax = p.Results.Nmax;
      end

      n=[1:Nmax];

      m = tmatrix.k_particle/tmatrix.k_medium;

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

      if p.Results.internal == false
        % Calculate external T-matrix
        b = -( j1d.*j0 - m*j0d.*j1 ) ./ ( j1d.*h0 - m*h0d.*j1 );
        a = -( j0d.*j1 - m*j1d.*j0 ) ./ ( h0d.*j1 - m*j1d.*h0 );
        T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
            [a(indexing);b(indexing)]);
      else
        % Calculate internal T-matrix
        d = -( h0d.*j0 - j0d.*h0 ) ./ ( m*j1.*h0d - j1d.*h0 );
        c = -( j0d.*h0 - h0d.*j0 ) ./ ( m*j1d.*h0 - h0d.*j1 );
        T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
            [c(indexing);d(indexing)]);
      end

      % Store T-matrix
      tmatrix.data = T;

      % TODO: We should store if it is internal or external
      % TODO: Is this a total or scattered T-matrix?  Store this too.
    end
  end
end
