classdef HermiteGaussian < ott.beam.BscFinite ...
    & ott.beam.properties.HermiteGaussian
% Construct a VSWF representation of tightly focussed Hermite-Gaussian beam.
% Inherits from :class:`ott.beam.BscFinite` and
% :class:`ott.beam.properties.HermiteGaussian`.
%
% Properties
%   - waist         -- Beam waist radius [m]
%   - index_medium  -- Refractive index of the medium
%   - omega         -- Optical angular frequency of light [1/s]
%   - position      -- Position of the beam [m]
%   - rotation      -- Rotation of the beam [3x3 rotation matrix]
%   - power         -- Beam power [W]
%   - mmode         -- Hermite mode number
%   - nmode         -- Hermite mode number
%   - polbasis      -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield      -- (2 numeric) Field in theta/phi or x/y directions
%   - mapping       -- Paraxial to far-field beam mapping
%   - data          -- Internal BSC instance describing beam
%
% Methods
%   - getData       -- Get data for specific Nmax
%   - recalculate   -- Recalculate the beam data
%
% Static methods
%   - FromNa        -- Construct a beam specifying NA instead of waist

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function beam = FromNa(Na, varargin)
      % Construct a Hermite-Gaussian beam specifying NA instead of beam waist.
      %
      % Usage
      %   beam = HermiteGaussian.FromNa(NA, ...)
      %
      % For other parameters, see class constructor.

      p = inputParser;
      p.addOptional('mmode', 0, @isnumeric);
      p.addOptional('nmode', 0, @isnumeric);
      p.addParameter('mapping', 'sin');
      p.addParameter('index_medium', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      paraxial_order = p.Results.mmode + p.Results.nmode;

      waist = ott.bsc.LgParaxial.WaistFromNa(Na, p.Results.index_medium, ...
          paraxial_order, p.Results.mapping);

      beam = ott.beam.HermiteGaussian(waist, varargin{:});
    end
  end

  methods
    function beam = HermiteGaussian(varargin)
      % Construct a VSWF Harmite-Gaussian beam representation.
      %
      % Usage
      %   beam = HermiteGaussian(waist, mmode, nmode, ...)
      %
      % Optional named arguments
      %   - waist (numeric) -- Beam waist [m].  Default: ``1e-6``.
      %
      %   - mmode (numeric) -- Mode number.  Default: ``0``.
      %
      %   - nmode (numeric) -- Mode number.  Default: ``0``.
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
      %
      %   - calculate (logical) -- If the beam should be calculated after
      %     construction.  Set to false if you intend to change beam
      %     properties immediately after construction.
      %     Uses :meth:`recalculate`.  Default: ``true``.

      p = inputParser;
      p.addOptional('waist', 1.0e-6, @isnumeric);
      p.addOptional('mmode', 0, @isnumeric);
      p.addOptional('nmode', 0, @isnumeric);
      p.addParameter('polfield', [1, 1i], @isnumeric);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('mapping', 'sin');
      p.addParameter('index_medium', 1.0);
      p.addParameter('power', 1.0);
      p.addParameter('calculate', true);
      p.parse(varargin{:});

      beam.waist = p.Results.waist;
      beam.mmode = p.Results.mmode;
      beam.nmode = p.Results.nmode;
      beam.polfield = p.Results.polfield;
      beam.polbasis = p.Results.polbasis;
      beam.mapping = p.Results.mapping;
      beam.index_medium = p.Results.index_medium;
      beam.power = p.Results.power;

      if p.Results.calculate
        beam = beam.recalculate([]);
      end
    end

    function beam = recalculate(beam, ~)
      % Re-calculate data for beam.
      %
      % This function can be called to pre-compute the beam data for
      % the current beam parameters.  It is automatically called when
      % the beam is used, however, calling the function explicitly can
      % speed up run-time.
      %
      % Usage
      %   beam = beam.recalcualte(Nmax)
      %
      % Parameters
      %   - Nmax (numeric) -- Parameter is ignored.

      beam.data = ott.bsc.LgParaxial.FromHgMode(beam.waist, ...
          beam.mmode, beam.nmode, beam.polbasis, beam.polfield);
      beam.data = beam.data ./ beam.data.power .* beam.power;
    end
  end
end

