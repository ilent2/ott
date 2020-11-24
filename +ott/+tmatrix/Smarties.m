classdef Smarties < ott.tmatrix.Homogeneous
% Constructs a T-matrix using SMARTIES.
% Inherits from :class:`ott.tmatrix.Homogeneous`.
%
% SMARTIES is a method for calculating T-matrices for spheroids, full
% details can be found in
%
%   W.R.C. Somerville, B. Auguié, E.C. Le Ru,
%   JQSRT, Volume 174, May 2016, Pages 39-55.
%   https://doi.org/10.1016/j.jqsrt.2016.01.005
%
% SMARTIES is distributed with a Creative Commons Attribution-NonCommercial
% 4.0 International License.  Copyright 2015 Walter Somerville,
% Baptiste AuguiÃ©, and Eric Le Ru.
% This version of OTT includes a minimal version of SMARTIES for T-matrix
% calculation, if you use the SMARTIES code, please cite the above paper.
%
% Properties
%   - ordinary         -- Ordinary radius
%   - extraordinary    -- Extra-ordinary radius
%   - index_relative   -- Relative refractive index of particle
%
% Static methods
%   - FromShape       -- Construct a T-matrix from a shape description
%
% See also Smarties, :class:`ott.tmatrix.Mie` and :mod:`ott.tmatrix.smarties`.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    ordinary          % Ordinary radius
    extraordinary     % Extra-ordinary radius
  end

  methods (Static)
    function varargout = FromShape(shape, varargin)
      % Construct a T-matrix using SMARTIES/EBCM for spheroids.
      %
      % Usage
      %   tmatrix = Smarties.FromShape(shape, index_relative, ...)
      %
      %   [texternal, tinternal] = Smarties.FromShape(...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- A spheroid with the extraordinary
      %     axis aligned to the z-axis.
      %
      %   - index_relative (numeric) -- The relative refractive index
      %     of the particle compared to the surrounding medium.
      %
      % All other parameters passed to class constructor.

      assert(numel(shape) == 1 && isa(shape, 'ott.shape.Shape'), ...
          'shape must be a single ott.shape.Shape');
      assert(shape.zRotSymmetry == 0, ...
          'shape must be rotationally symmetric about z-axis');

      if isa(shape, 'ott.shape.Sphere')
        [varargout{1:nargout}] = ott.tmatrix.Smarties(shape.radius, ...
            shape.radius, varargin{:});
      elseif isa(shape, 'ott.shape.Ellipsoid')
        [varargout{1:nargout}] = ott.tmatrix.Smarties(shape.radii(1), ...
            shape.radii(3), varargin{:});
      elseif isa(shape, 'ott.shape.Superellipsoid')
        assert(shape.isEllipsoid, 'Shape must be an spheroid');
        [varargout{1:nargout}] = ott.tmatrix.Smarties(shape.radii(1), ...
            shape.radii(3), varargin{:});
      else
        error('Shape must be a sphere, ellipsoid or supeperellipsoid');
      end
    end
  end

  methods
    function [tmatrix, iTmatrix] = Smarties(varargin)
      % Construct a T-matrix using SMARTIES/EBCM for spheroids.
      %
      % Usage
      %   tmatrix = Smarties(ordinary, extraordinary, index_relative...)
      %
      %   [texternal, tinternal] = Smarties(...)
      %
      % Parameters
      %   - ordinary (numeric) -- Ordinary radius.
      %
      %   - extraordinary (numeric) -- Extraordinary radius.
      %
      %   - index_relative (numeric) -- The relative refractive index
      %     of the particle compared to the surrounding medium.
      %
      % Optional named parameters
      %   - internal (logical) -- If true, the returned T-matrix is
      %     an internal T-matrix.  Ignored for two outputs.
      %     Default: ``false``.
      %
      %   - Nmax (numeric) -- Size of the VSWF expansion used for the
      %     T-matrix calculation.
      %     Default: ``ott.utis.ka2nmax(2*pi*max(radius))`` (external) or
      %     ``ott.utis.ka2nmax(2*pi*max(radius)*index_relative)`` (internal).
      %
      %   - npts (numeric) -- Number of points for surface integral.
      %     Default: ``Nmax*Nmax``.
      %
      %   - verbose (logical) -- Enables additional output from SMARTIES.
      %     Default: ``false``.

      p = inputParser;
      p.addOptional('ordinary', @isnumeric);
      p.addOptional('extraordinary', @isnumeric);
      p.addOptional('index_relative', 1.0, @isnumeric);
      p.addParameter('internal', false);
      p.addParameter('Nmax', []);
      p.addParameter('npts', []);
      p.addParameter('verbose', false);
      p.parse(varargin{:});

      % Store parameters
      tmatrix.ordinary = p.Results.ordinary;
      tmatrix.extraordinary = p.Results.extraordinary;
      tmatrix.index_relative = p.Results.index_relative;

      % Get/Check Nmax
      Nmax = ott.tmatrix.Tmatrix.getValidateNmax(...
          p.Results.Nmax, max(tmatrix.ordinary, tmatrix.extraordinary), ...
          tmatrix.index_relative, p.Results.internal || nargout ~= 1);

      % Get or estimate npts from inputs
      npts = p.Results.npts;
      if isempty(npts)
        npts = Nmax.^2;
      else
        assert(isnumeric(npts) && isscalar(npts) && npts > 0, ...
            'npts must be positive numeric scalar');
      end

      % Setup structure for SMARTIES inputs
      stParams = struct();
      stParams.a = tmatrix.ordinary;
      stParams.c = tmatrix.extraordinary;
      stParams.k1 = 2*pi;
      stParams.s = tmatrix.index_relative;
      stParams.N = Nmax;
      stParams.nNbTheta = npts;

      % Setup structure of optional parameters
      stOptions.bGetR = nargout == 2 || p.Results.internal;
      stOptions.Delta = 0;
      stOptions.NB = 0;
      stOptions.bGetSymmetricT = false;
      stOptions.bOutput = p.Results.verbose;

      % Solve for the T-matrix
      [~, CstTRa] = ott.tmatrix.smarties.slvForT(stParams, stOptions);
      Texternal = tmatrix.getTmatrixData(CstTRa, 'st4MT');
      if nargout == 2 || p.Results.internal
        Tinternal = tmatrix.getTmatrixData(CstTRa, 'st4MR');
      end

      % Assign outputs
      if nargout == 2

        iTmatrix = tmatrix;
        iTmatrix.data = Tinternal;
        iTmatrix = iTmatrix.setType('internal');

        tmatrix.data = Texternal;
        tmatrix = tmatrix.setType('scattered');

      else

        if p.Results.internal
          tmatrix.data = Tinternal;
          tmatrix = tmatrix.setType('internal');
        else
          tmatrix.data = Texternal;
          tmatrix = tmatrix.setType('scattered');
        end

      end
    end
  end

  methods (Hidden, Static)
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

  methods % Getters/setters
    function tmatrix = set.ordinary(tmatrix, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'ordinary radius must be positive numeric scalar');
      tmatrix.ordinary = val;
    end

    function tmatrix = set.extraordinary(tmatrix, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'extraordinary radius must be positive numeric scalar');
      tmatrix.extraordinary = val;
    end
  end
end
