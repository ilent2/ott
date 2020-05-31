classdef InceGaussian < ott.beam.vswf.BscScalar ...
    & ott.beam.properties.InceGaussian ...
    & ott.beam.properties.FarfieldMapping
% Ince-Gaussian beam VSWF representation using LG-paraxial point matching.
% Inherits from :class:`ott.beam.vswf.BscScalar` and
% :class:`ott.beam.abstract.InceGaussian`.
%
% Properties
%   - waist         -- Beam waist radius
%   - lmode         -- Azimuthal Laguerre mode order
%   - porder        -- Radial Laguerre mode order
%   - ellipticity   -- Ellipticity of beam
%   - parity        -- Parity of beam
%   - medium        -- Properties of the optical medium
%   - omega         -- Optical frequency of light
%   - position      -- Position of the beam
%   - rotation      -- Rotation of the beam
%   - power         -- Beam power
%   - polarisation  -- Far-field polarisation of paraxial beam
%   - mapping       -- Paraxial to far-field mapping

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.properties.FarfieldMapping.likeProperties(other, args);
      args = ott.beam.vswf.BscScalar.likeProperties(other, args);
      args = ott.beam.properties.InceGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Construct a VSWF beam like another beam
      args = ott.beam.vswf.InceGaussian.likeProperties(other, varargin);
      beam = ott.beam.vswf.InceGaussian(args{:});
    end
  end

  methods
    function bsc = InceGaussian(varargin)
      % Construct a VSWF representation of a Ince-Gaussian beam
      %
      % Usage
      %   bsc = InceGaussian(waist, lmode, porder, parity, ellipticity, ...)
      %
      % Parameters
      %   - waist (numeric) -- Beam waist.
      %   - lmode (numeric) -- Azimuthal mode number.
      %   - porder (numeric) -- Paraxial mode order.
      %   - parity (enum) -- Either 'even' or 'odd'.
      %   - ellipticity (numeric) -- Ellipticity of coordinates.
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
      %     Can be either 'sin' or 'tan' (small angle).
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).  Default: ``'sin'``.

      p = inputParser;
      p.addOptional('waist', [], @isnumeric);
      p.addOptional('lmode', [], @isnumeric);
      p.addOptional('porder', [], @isnumeric);
      p.addOptional('parity', [], @(x) any(strcmpi(x, {'odd', 'even'})));
      p.addOptional('ellipticity', [], @isnumeric);
      p.addParameter('mapping', 'sin');
      p.addParameter('polarisation', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bsc = bsc@ott.beam.properties.FarfieldMapping(...
          p.Results.mapping, 'hemisphere', 'neg');
      bsc = bsc@ott.beam.properties.InceGaussian(...
          p.Results.waist, p.Results.lmode, p.Results.porder, ...
          p.Results.parity, p.Results.ellipticity, ...
          'polarisation', p.Results.polarisation, unmatched{:});
      bsc = bsc@ott.beam.vswf.BscScalar();

      % Construct a beam LgParaxialBasis beam
      data = ott.beam.vswf.LgParaxialBasis.FromIgMode(...
          bsc.waist, bsc.lmode, bsc.porder, ...
          bsc.parity, bsc.ellipticity, ...
          bsc.polarisation, 'mapping', bsc.mapping);
      bsc = bsc.setCoefficients(sum(data));
    end
  end
end
