classdef LaguerreGaussian < ott.beam.BscFinite ...
    & ott.beam.properties.LaguerreGaussian
% Construct a VSWF representation of tightly focussed Laguerre-Gaussian beam.
% Inherits from :class:`ott.beam.BscFinite` and
% :class:`ott.beam.properties.LaguerreGaussian`.
%
% Properties
%   - waist         -- Beam waist radius [m]
%   - index_medium  -- Refractive index of the medium
%   - omega         -- Optical angular frequency of light [1/s]
%   - position      -- Position of the beam [m]
%   - rotation      -- Rotation of the beam [3x3 rotation matrix]
%   - power         -- Beam power [W]
%   - lmode         -- Azimuthal mode number
%   - pmode         -- Radial mode number
%   - polbasis      -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- (2 numeric) Field in theta/phi or x/y directions
%   - mapping       -- Paraxial to far-field beam mapping
%   - data          -- Internal BSC instance describing beam
%
% Methods
%   - getData       -- Get data for specific Nmax
%   - recalculate   -- Recalculate the beam data

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = LaguerreGaussian(varargin)
      % Construct a VSWF Laguerre-Gaussian beam represntation.
      %
      % Usage
      %   beam = LaguerreGaussian(waist, lmode, pmode, ...)
      %
      % Optional named arguments
      %   - waist (numeric) -- Beam waist [m].  Default: ``1e-6``.
      %
      %   - lmode (numeric) -- Azimuthal mode number.  Default: ``0``.
      %
      %   - pmode (numeric) -- Radial mode number.  Default: ``0``.
      %
      %   - polbasis (enum) -- Polarisation basis.  Either 'polar' or
      %     'cartesian'.  Default: ``'cartesian'``.
      %
      %   - polfield (2 numeric) -- Field in either the theta/phi or
      %     x/y directions (depending on basis).  Default: ``[1, 1i]``.
      %
      %   - mapping (enum) -- Mapping method for paraxial far-field.
      %     Can be either 'sin', 'tan' (small angle) or 'theta'.
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).  Default: ``'sin'``.
      %
      %   - index_medium (numeric) -- Refractive index of the medium.
      %     Default: ``1.0``.
      %
      %   - power (numeric) -- Beam power [W].  Default: ``1.0``.

      p = inputParser;
      p.addOptional('waist', 1.0e-6, @isnumeric);
      p.addOptional('lmode', 0, @isnumeric);
      p.addOptional('pmode', 0, @isnumeric);
      p.addParameter('polfield', [1, 1i], @isnumeric);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('mapping', 'sin');
      p.addParameter('index_medium', 1.0);
      p.addParameter('power', 1.0);
      p.parse(varargin{:});

      beam.waist = p.Results.waist;
      beam.lmode = p.Results.lmode;
      beam.pmode = p.Results.pmode;
      beam.polfield = p.Results.polfield;
      beam.polbasis = p.Results.polbasis;
      beam.mapping = p.Results.mapping;
      beam.index_medium = p.Results.index_medium;
      beam.power = p.Results.power;

      beam = beam.recalculate([]);
    end
  end

  methods (Hidden)
    function [data, vswfData] = recalculateInternal(beam, Nmax, vswfData)
      % Re-calculate BSC data for specified Nmax.

      % TODO
    end
  end
end

