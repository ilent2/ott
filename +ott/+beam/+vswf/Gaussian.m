classdef Gaussian < ott.beam.vswf.BscScalar ...
    & ott.beam.properties.Gaussian ...
    & ott.beam.utils.FarfieldMapping
% Gaussian beam VSWF representation using LG-paraxial point matching.
% Inherits from :class:`ott.beam.vswf.BscScalar`,
% :class:`ott.beam.properties.Gaussian`.
%
% Properties
%   - waist         -- Beam waist radius
%   - medium        -- Properties of the optical medium
%   - omega         -- Optical frequency of light
%   - position      -- Position of the beam
%   - rotation      -- Rotation of the beam
%   - power         -- Beam power
%   - polarisation  -- Far-field polarisation of paraxial beam
%   - mapping       -- Paraxial to Far-field mapping

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.FarfieldMapping.likeProperties(other, args);
      args = ott.beam.vswf.BscScalar.likeProperties(other, args);
      args = ott.beam.properties.Gaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Construct a VSWF beam like another beam
      args = ott.beam.vswf.Gaussian.likeProperties(other, varargin);
      beam = ott.beam.vswf.Gaussian(args{:});
    end
  end

  methods
    function bsc = Gaussian(varargin)
      % Construct a VSWF representation of a Gaussian beam
      %
      % Usage
      %   bsc = Gaussian(waist, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist.
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
      p.addParameter('mapping', 'sintheta');
      p.addParameter('polarisation', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bsc = bsc@ott.beam.utils.FarfieldMapping(p.Results.mapping);
      bsc = bsc@ott.beam.properties.Gaussian(p.Results.waist, ...
          'polarisation', p.Results.polarisation, unmatched{:});
      bsc = bsc@ott.beam.vswf.BscScalar();

      % Construct a beam LgParaxialBasis beam
      % TODO: Perhaps this should be called later?  Grow?
      data = ott.beam.vswf.LgParaxialBasis.FromLgMode(...
          bsc.waist, 0, 0, bsc.polarisation, 'mapping', bsc.mapping);
      bsc = bsc.setCoefficients(sum(data));
    end
  end
end
