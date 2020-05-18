classdef LaguerreGaussian < ott.beam.vswf.BscScalar ...
    & ott.beam.abstract.LaguerreGaussian ...
% Laguerre-Gaussian beam VSWF representation using LG-paraxial point matching.
% Inherits from :class:`ott.beam.vswf.BscScalar` and
% :class:`ott.beam.abstract.LaguerreGaussian`.
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

  properties
    polarisation      % Far-field polarisation of paraxial beam
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

      p = inputParser;
      p.addOptional('waist', 1.0, @isnumeric);
      p.addOptional('lmode', 0, @isnumeric);
      p.addOptional('pmode', 0, @isnumeric);
      p.addParameter('polarisation', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct a beam LgParaxialBasis beam
      data = ott.beam.vswf.LgParaxialBasis.FromLgMode(p.Results.waist, ...
          p.Results.lmode, p.Results.pmode, p.Results.polarisation);

      % Construct this object from data
      bsc = bsc@ott.beam.vswf.BscScalar(sum(data));
      bsc = bsc@ott.beam.abstract.LaguerreGaussian(p.Results.waist, ...
          p.Results.lmode, p.Results.pmode, unmatched{:});
      
      bsc.polarisation = p.Results.polarisation;
    end
  end
  
  methods (Hidden)
    function bsc = setBeamPower(bsc, power)
      bsc = setBeamPower@ott.beam.vswf.Bsc(bsc, power);
      bsc = setBeamPower@ott.beam.abstract.LaguerreGaussian(bsc, power);
    end
    
    function p = getBeamPower(bsc)
      p = getBeamPower@ott.beam.vswf.Bsc(bsc);
    end
  end
end
