classdef Gaussian < ott.beam.BscFinite ...
    & ott.beam.properties.Gaussian
% Construct a VSWF representation of a tightly focussed Gaussian beam.
% Inherits from :class:`ott.beam.BscFinite` and
% :class:`ott.beam.properties.Gaussian`.
%
% Properties
%   - waist         -- Beam waist radius [m]
%   - index_medium  -- Refractive index of the medium
%   - omega         -- Optical angular frequency of light [1/s]
%   - position      -- Position of the beam [m]
%   - rotation      -- Rotation of the beam [3x3 rotation matrix]
%   - power         -- Beam power [W]
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
    function beam = Gaussian(varargin)
      % Construct a VSWF representation of a Gaussian beam
      %
      % Usage
      %   beam = Gaussian(waist, ...)
      %
      % Optional named arguments
      %   - waist (numeric) -- Beam waist [m].  Default: ``1e-6``.
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
      p.addOptional('waist', 1e-6);
      p.addParameter('polfield', [1, 1i]);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('mapping', 'sin');
      p.addParameter('index_medium', 1.0);
      p.addParameter('power', 1.0);
      p.parse(varargin{:});

      beam = beam.recalculate([]);
    end
  end
  
  methods (Hidden)
    function [data, vswfData] = recalculateInternal(beam, ~, vswfData)
      % Re-calculate BSC data for specified Nmax.

      % TODO
    end
  end
end
