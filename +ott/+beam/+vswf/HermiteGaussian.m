classdef HermiteGaussian < ott.beam.vswf.BscScalar ...
    & ott.beam.properties.HermiteGaussian ...
    & ott.beam.utils.FarfieldMapping
% Hermite-Gaussian beam VSWF representation using LG-paraxial point matching.
% Inherits from :class:`ott.beam.vswf.BscScalar`,
% :class:`ott.beam.abstract.HermiteGaussian`.
%
% Properties
%   - waist         -- Beam waist radius
%   - mmode         -- Azimuthal Laguerre mode order
%   - nmode         -- Radial Laguerre mode order
%   - medium        -- Properties of the optical medium
%   - omega         -- Optical frequency of light
%   - position      -- Position of the beam
%   - rotation      -- Rotation of the beam
%   - power         -- Beam power
%   - mapping       -- Paraxial to Far-field mapping

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.FarfieldMapping.likeProperties(other, args);
      args = ott.beam.vswf.BscScalar.likeProperties(other, args);
      args = ott.beam.properties.HermiteGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Construct a VSWF beam like another beam
      args = ott.beam.vswf.HermiteGaussian.likeProperties(other, varargin);
      beam = ott.beam.vswf.HermiteGaussian(args{:});
    end
  end

  methods
    function bsc = HermiteGaussian(varargin)
      % Construct a VSWF representation of a Laguerre-Gaussian beam
      %
      % Usage
      %   bsc = LaguerreGaussian(waist, mmode, nmode, ...)
      %
      %   bsc = LaguerreGaussian(...)
      %   Construct a Gaussian beam with default parameters.
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist.  Default: ``1.0``.
      %   - mmode (numeric) -- Mode number.  Default: ``0``.
      %   - nmode (numeric) -- Mode number.  Default: ``0``.
      %
      % Optional named parameters
      %   - position (3x1 numeric) -- Position of the beam.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Rotation of the beam.
      %     Default: ``eye(3)``.
      %
      %   - polarisation (2-numeric) -- Polarisation of the beam.
      %     2 element Jones vector.  Default: ``[1, 0]``.
      %
      %   - mapping (enum) -- Mapping method for paraxial far-field.
      %     Can be either 'sintheta' or 'tantheta' (small angle).
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).  Default: ``'sintheta'``.

      p = inputParser;
      p.addOptional('waist', [], @isnumeric);
      p.addOptional('mmode', [], @isnumeric);
      p.addOptional('nmode', [], @isnumeric);
      p.addParameter('mapping', 'sintheta');
      p.addParameter('polarisation', [1, 0]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bsc = bsc@ott.beam.utils.FarfieldMapping(p.Results.mapping);
      bsc = bsc@ott.beam.properties.HermiteGaussian(...
          p.Results.waist, p.Results.mmode, p.Results.nmode, ...
          'polarisation', p.Results.polarisation, unmatched{:});
      bsc = bsc@ott.beam.vswf.BscScalar();

      % Construct a beam LgParaxialBasis beam
      data = ott.beam.vswf.LgParaxialBasis.FromHgMode(...
          bsc.waist, bsc.mmode, bsc.nmode, ...
          bsc.polarisation, 'mapping', bsc.mapping);
      bsc = bsc.setCoefficients(sum(data));
    end
  end
end
