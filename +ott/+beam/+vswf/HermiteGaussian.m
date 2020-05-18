classdef HermiteGaussian < ott.beam.vswf.BscScalar ...
    & ott.beam.abstract.HermiteGaussian ...
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

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    polarisation      % Far-field polarisation of paraxial beam
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
      %     2 element Jones vector.  Default: ``[1, 1i]``.

      p = inputParser;
      p.addOptional('waist', 1.0, @isnumeric);
      p.addOptional('mmode', 0, @isnumeric);
      p.addOptional('nmode', 0, @isnumeric);
      p.addParameter('polarisation', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct a beam LgParaxialBasis beam
      data = ott.beam.vswf.LgParaxialBasis.FromHgMode(...
          p.Results.waist, p.Results.mmode, p.Results.nmode, ...
          p.Results.polarisation);

      % Construct this object from data
      bsc = bsc@ott.beam.vswf.BscScalar(sum(data));
      bsc = bsc@ott.beam.abstract.HermiteGaussian(p.Results.waist, ...
          p.Results.mmode, p.Results.nmode, unmatched{:});
      
      bsc.polarisation = p.Results.polarisation;
    end
  end
  
  methods (Hidden)
    function bsc = setBeamPower(bsc, power)
      bsc = setBeamPower@ott.beam.vswf.Bsc(bsc, power);
      bsc = setBeamPower@ott.beam.abstract.HermiteGaussian(bsc, power);
    end
    
    function p = getBeamPower(bsc)
      p = getBeamPower@ott.beam.vswf.Bsc(bsc);
    end
  end
end
