classdef LaguerreGaussian < ott.beam.vswf.BscScalar ...
    & ott.beam.properties.LaguerreGaussian ...
    & ott.beam.utils.FarfieldMapping
% Laguerre-Gaussian beam VSWF representation using LG-paraxial point matching.
% Inherits from :class:`ott.beam.vswf.BscScalar` and
% :class:`ott.beam.properties.LaguerreGaussian`.
%
% This is a simplified interface for :class:`LgParaxialBasis`.
%
% Properties
%   - waist         -- Beam waist radius
%   - lmode         -- Azimuthal Laguerre mode order
%   - pmode         -- Radial Laguerre mode order
%   - medium        -- Properties of the optical medium
%   - omega         -- Optical frequency of light
%   - position      -- Position of the beam
%   - rotation      -- Rotation of the beam
%   - power         -- Beam power
%   - polarisation  -- Far-field polarisation of paraxial beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.FarfieldMapping.likeProperties(other, args);
      args = ott.beam.vswf.BscScalar.likeProperties(other, args);
      args = ott.beam.properties.LaguerreGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Construct a VSWF beam like another beam
      args = ott.beam.vswf.LaguerreGaussian.likeProperties(other, varargin);
      beam = ott.beam.vswf.LaguerreGaussian(args{:});
    end
  end

  methods
    function bsc = LaguerreGaussian(varargin)
      % Construct a VSWF representation of a Laguerre-Gaussian beam
      %
      % Usage
      %   bsc = LaguerreGaussian(waist, lmode, pmode, ...)
      %
      %   bsc = LaguerreGaussian(...)
      %   Construct a Gaussian beam with default parameters.
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist.  Default: ``1.0``.
      %   - lmode (numeric) -- Azimuthal mode.  Default: ``0``.
      %   - pmode (numeric) -- Radial mode. Default: ``0``.
      %
      % Optional named parameters
      %   - position (3x1 numeric) -- Position of the beam.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Rotation of the beam.
      %     Default: ``eye(3)``.
      %
      %   - polarisation (2-numeric) -- Polarisation of the beam.
      %     2 element Jones vector.  Default: ``[1, 1i]``.
      %
      %   - mapping (enum) -- Mapping method for paraxial far-field.
      %     Can be either 'sintheta' or 'tantheta' (small angle).
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).  Default: ``'sintheta'``.

      p = inputParser;
      p.addOptional('waist', [], @isnumeric);
      p.addOptional('lmode', [], @isnumeric);
      p.addOptional('pmode', [], @isnumeric);
      p.addParameter('mapping', 'sintheta');
      p.addParameter('polarisation', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bsc = bsc@ott.beam.utils.FarfieldMapping(p.Results.mapping);
      bsc = bsc@ott.beam.properties.LaguerreGaussian(...
          p.Results.waist, p.Results.lmode, p.Results.pmode, ...
          'polarisation', p.Results.polarisation, unmatched{:});
      bsc = bsc@ott.beam.vswf.BscScalar();

      % Construct a beam LgParaxialBasis beam
      data = ott.beam.vswf.LgParaxialBasis.FromLgMode(bsc.waist, ...
          bsc.lmode, bsc.pmode, bsc.polarisation, ...
          'mapping', bsc.mapping);
      bsc = bsc.setCoefficients(sum(data));
    end
  end
end
