classdef TmatrixSmarties < ott.Tmatrix
% Constructs a T-matrix using SMARTIES.
% Inherits from :class:`+ott.Tmatrix`
%
% SMARTIES is a method for calculating T-matrices for spheroids.
%
% This class requires the SMARTIES toolbox has been loaded onto the path.
% Details about the toolbox can be found in:
%
%   Somerville, Auguié, Le Ru.  JQSRT, Volume 174, May 2016, Pages 39-55.
%   https://doi.org/10.1016/j.jqsrt.2016.01.005
%
% See also TmatrixSmarties

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties %(SetAccess=protected)
    k_medium            % Wavenumber of medium
    k_particle          % Wavenumber of particle
  end

  methods (Static)
    function tmatrix = simple(shape, varargin)
      % Construct a T-matrix using SMARTIES/EBCM for spheroids
      %
      % simple(shape) constructs a new simple T-matrix for the given shape.
      %
      % simple(name, parameters) constructs a new T-matrix for the
      % shape described by name and parameters.
      %
      % Shape must be spheroid with extraordinary axis aligned to z-axis.
      %
      % Optional named parameters:
      %     'internal'    bool    Calculate internal T-matrix (default: 0)

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('parameters', []);
      p.parse(varargin{:});

      % Get a shape object from the inputs
      if ischar(shape) && ~isempty(p.Results.parameters)
        shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
        varargin = varargin(2:end);
      elseif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
        error('Must input either Shape object or string and parameters');
      end

      % Check that the particle is rotationally symmetric
      axsym = shape.axialSymmetry();
      if axsym(3) ~= 0
        error('Only axially symetric particles supported for now');
      end

      % Construct the T-matrix
      if isa(shape, 'ott.shapes.Sphere')
        tmatrix = ott.TmatrixSmarties(shape.radius, ...
            shape.radius, varargin{:});
      elseif isa(shape, 'ott.shapes.Ellipsoid')
        tmatrix = ott.TmatrixSmarties(shape.a, ...
            shape.c, varargin{:});
      else
        error('Only supports sphere and ellipsoid shapes for now');
      end
    end
  end

  methods
    function tmatrix = TmatrixSmarties(a, c, varargin)
      % Construct a T-matrix using SMARTIES/EBCM for spheroids
      %
      % TmatrixSmarties(a, c, ...) constructs the T-matrix for a
      % spheroid with ordinary radius a and extraordinary radius c.
      %
      % Optional named parameters:
      %     internal             bool    Calculate internal T-matrix
      %         Default: false
      %
      %     k_medium             num     Wavenumber in medium
      %     wavelength_medium    num     Wavelength in medium
      %     index_medium         num     Refractive index in medium
      %         Default: k_medium = 2*pi
      %
      %     k_particle           num     Wavenumber in particle
      %     wavelength_particle  num     Wavelength in particle
      %     index_particle       num     Refractive index in particle
      %         Default: k_particle = 2*pi*index_relative
      %
      %     index_relative       num     Relative refractive index
      %         Default: 1.0
      %     wavelength0          num     Wavelength in vacuum
      %         Default: 1.0
      %
      %     Nmax                 int     Maximum multiple order for
      %         (default: ott.utils.ka2nmax(max_radius*wavenumber))
      %     npts                 int     Number of points to use for
      %         surface integral (default: Nmax*Nmax)
      %     verbose              bool    Display additional output

      tmatrix = tmatrix@ott.Tmatrix();

      % Parse inputs
      p = inputParser;

      p.addParameter('internal', false);

      p.addParameter('wavelength0', []);
      p.addParameter('index_relative', []);

      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);

      p.addParameter('k_particle', []);
      p.addParameter('wavelength_particle', []);
      p.addParameter('index_particle', []);

      p.addParameter('Nmax', []);
      p.addParameter('npts', []);
      p.addParameter('verbose', false);
      p.addParameter('progress_callback', @(x) []);

      p.parse(varargin{:});
      
      % First check we have SMARTIES
      if exist('slvForT', 'file') ~= 2
        error('Could not find SMARTIES function slvForT, check SMARTIES installed');
      end

      % Get the wavenumber of the particle and medium
      [tmatrix.k_medium, tmatrix.k_particle] = ...
          tmatrix.parser_wavenumber(p, 2*pi);

      % Get or estimate Nmax from the inputs
      Nmax = p.Results.Nmax;
      if isempty(Nmax)
        maxRadius = max(a, c);
        if p.Results.internal
          Nmax = ott.utils.ka2nmax(maxRadius * abs(tmatrix.k_particle));
        else
          Nmax = ott.utils.ka2nmax(maxRadius * abs(tmatrix.k_medium));
        end
      end

      % Get or estimate npts from inputs
      npts = p.Results.npts;
      if isempty(npts)
        npts = Nmax.^2;
      end

      % Setup structure of simulation inputs
      stParams.a = a .* tmatrix.k_medium;
      stParams.c = c .* tmatrix.k_medium;
      stParams.k1 = 1.0;
      stParams.s = tmatrix.k_particle./tmatrix.k_medium;
      stParams.N = Nmax;
      stParams.nNbTheta = npts;

      % Setup structure of optional parameters
      stOptions.bGetR = p.Results.internal;
      stOptions.Delta = 0;
      stOptions.NB = 0;
      stOptions.bGetSymmetricT = false;
      stOptions.bOutput = p.Results.verbose;

      % Solve for the T-matrix
      [~, CstTRa] = slvForT(stParams, stOptions);

      % Store the result
      if p.Results.internal
        tmatrix.type = 'internal';
        tmatrix.data = ott.TmatrixSmarties.getTmatrixData(CstTRa, 'st4MR');
      else
        tmatrix.type = 'scattered';
        tmatrix.data = ott.TmatrixSmarties.getTmatrixData(CstTRa, 'st4MT');
      end
    end
  end

  methods (Access=protected, Static)
    function data = getTmatrixData(CstTRa, name)
      % Convert the SMARTIES structure to T-matrix data

      Nmax = length(CstTRa)-1;
      tsz = ott.utils.combined_index(Nmax, Nmax);
      data = sparse(2*tsz, 2*tsz);
      sz = size(data);

      % Loop over each multipole solution
      for m = -Nmax:Nmax

        n = abs(m)-1;
        if m == 0, n = 0; end

        % First do oe
        blk = CstTRa{abs(m)+1}.([name, 'oe']);
        assert(abs(m) == blk.m);

        ridx = ott.utils.combined_index(n+blk.ind1, m);
        cidx = ott.utils.combined_index(n+blk.ind1, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx, cidx)) = blk.M11;

        ridx = ott.utils.combined_index(n+blk.ind1, m);
        cidx = ott.utils.combined_index(n+blk.ind2, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx, cidx+tsz)) = blk.M12;

        ridx = ott.utils.combined_index(n+blk.ind2, m);
        cidx = ott.utils.combined_index(n+blk.ind1, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx+tsz, cidx)) = blk.M21;

        ridx = ott.utils.combined_index(n+blk.ind2, m);
        cidx = ott.utils.combined_index(n+blk.ind2, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx+tsz, cidx+tsz)) = blk.M22;

        % Now do eo
        blk = CstTRa{abs(m)+1}.([name, 'eo']);
        assert(abs(m) == blk.m);

        ridx = ott.utils.combined_index(n+blk.ind1, m);
        cidx = ott.utils.combined_index(n+blk.ind1, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx, cidx)) = blk.M11;

        ridx = ott.utils.combined_index(n+blk.ind1, m);
        cidx = ott.utils.combined_index(n+blk.ind2, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx, cidx+tsz)) = blk.M12;

        ridx = ott.utils.combined_index(n+blk.ind2, m);
        cidx = ott.utils.combined_index(n+blk.ind1, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx+tsz, cidx)) = blk.M21;

        ridx = ott.utils.combined_index(n+blk.ind2, m);
        cidx = ott.utils.combined_index(n+blk.ind2, m);
        [ridx, cidx] = meshgrid(ridx, cidx);
        data(sub2ind(sz, ridx+tsz, cidx+tsz)) = blk.M22;
      end
    end
  end
end
