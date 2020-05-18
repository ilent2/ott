classdef InceGaussian < ott.beam.vswf.BscScalar ...
    & ott.beam.abstract.InceGaussian ...
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

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    polarisation      % Far-field polarisation of paraxial beam
  end

  methods
    function bsc = InceGaussian(varargin)
      % Construct a VSWF representation of a Ince-Gaussian beam
      %
      % Usage
      %   bsc = InceGaussian(waist, lmode, porder, parity, ellipticity, ...)
      %
      %   bsc = InceGaussian(...)
      %   Construct a beam with `waist = 1`, `lmode = 0`, `porder = 1`,
      %   even parity and `ellipticity = 1`.
      %   Parameters can also be passed as named arguments.
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

      p = inputParser;
      p.addOptional('waist', 1.0, @isnumeric);
      p.addOptional('lmode', 0, @isnumeric);
      p.addOptional('porder', 1, @isnumeric);
      p.addOptional('parity', 'even', @(x) any(strcmpi(x, {'even', 'odd'})));
      p.addOptional('ellipticity', 1, @isnumeric);
      p.addParameter('polarisation', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct a beam LgParaxialBasis beam
      data = ott.beam.vswf.LgParaxialBasis.FromIgMode(...
          p.Results.waist, p.Results.lmode, p.Results.porder, ...
          p.Results.parity, p.Results.ellipticity, ...
          p.Results.polarisation);

      % Construct this object from data
      bsc = bsc@ott.beam.vswf.BscScalar(sum(data));
      bsc = bsc@ott.beam.abstract.InceGaussian(...
          p.Results.waist, p.Results.lmode, p.Results.porder, ...
          p.Results.parity, p.Results.ellipticity, unmatched{:});
      
      bsc.polarisation = p.Results.polarisation;
    end
  end
  
  methods (Hidden)
    function bsc = setBeamPower(bsc, power)
      bsc = setBeamPower@ott.beam.vswf.Bsc(bsc, power);
      bsc = setBeamPower@ott.beam.abstract.InceGaussian(bsc, power);
    end
    
    function p = getBeamPower(bsc)
      p = getBeamPower@ott.beam.vswf.Bsc(bsc);
    end
  end
end
