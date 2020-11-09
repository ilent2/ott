classdef InceGaussian < ott.beam.BscFinite ...
    & ott.beam.properties.InceGaussian
% Construct a VSWF representation of tightly focussed Laguerre-Gaussian beam.
% Inherits from :class:`ott.beam.BscFinite` and
% :class:`ott.beam.properties.InceGaussian`.
%
% Properties
%   - waist         -- Beam waist radius [m]
%   - index_medium  -- Refractive index of the medium
%   - omega         -- Optical angular frequency of light [1/s]
%   - position      -- Position of the beam [m]
%   - rotation      -- Rotation of the beam [3x3 rotation matrix]
%   - power         -- Beam power [W]
%   - lmode         -- Azimuthal mode number
%   - porder        -- Paraxial mode order.
%   - ellipticity   -- Ellipticity of coordinates
%   - parity        -- Parity of beam ('even' or 'odd')
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
      % Construct a InceGaussian beam specifying NA instead of beam waist.
      %
      % Usage
      %   beam = InceGaussian.FromNa(NA, ...)
      %
      % Optional named parameters
      %   - truncation_angle (numeric) -- Angle to truncate the beam.
      %     Defaults to match the beam NA: ``asin(NA/index_medium)``.
      %
      % For other parameters, see class constructor.

      p = inputParser;
      p.addOptional('lmode', 1, @isnumeric);    % Not used
      p.addOptional('porder', 0, @isnumeric);
      p.addParameter('mapping', 'sin');
      p.addParameter('index_medium', 1.0);
      p.addParameter('truncation_angle', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      waist = ott.bsc.LgParaxial.WaistFromNa(Na, p.Results.index_medium, ...
          p.Results.porder, p.Results.mapping);

      truncation_angle = p.Results.truncation_angle;
      if isempty(truncation_angle)
        truncation_angle = ott.utils.na2angle(Na, p.Results.index_medium);
      end

      beam = ott.beam.InceGaussian(waist, ...
          'truncation_angle', truncation_angle, varargin{:});
    end
  end

  methods
    function beam = InceGaussian(varargin)
      % Construct a VSWF Ince--Gaussian beam representation.
      %
      % Usage
      %   beam = InceGaussian(waist, lmode, porder, ellipticity, parity, ...)
      %
      % Optional named arguments
      %   - waist (numeric) -- Beam waist [m].  Default: ``1e-6``.
      %
      %   - lmode (numeric) -- Azimuthal mode order.  Default: ``1``.
      %
      %   - porder (numeric) -- Paraxial mode order.  Default: ``1``.
      %
      %   - ellipticity (numeric) -- Ellipticity of the beam.  Default: ``1``.
      %
      %   - parity (enum) -- Parity.  Either 'even' or 'odd'.
      %     Default: ``'even'``.
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
      %   - truncation_angle (numeric) -- Adds a hard edge truncation
      %     to the beam in the far-field.  Default: ``pi`` (i.e., no
      %     truncation added).
      %
      %   - calculate (logical) -- If the beam should be calculated after
      %     construction.  Set to false if you intend to change beam
      %     properties immediately after construction.
      %     Uses :meth:`recalculate`.  Default: ``true``.

      p = inputParser;
      p.addOptional('waist', 1.0e-6, @isnumeric);
      p.addOptional('lmode', 1, @isnumeric);
      p.addOptional('porder', 1, @isnumeric);
      p.addOptional('ellipticity', 1, @isnumeric);
      p.addOptional('parity', 'even', ...
          @(x) sum(strcmpi(x, {'even', 'odd'})) == 1);
      p.addParameter('polfield', [1, 1i], @isnumeric);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('mapping', 'sin');
      p.addParameter('power', 1.0);
      p.addParameter('truncation_angle', pi);
      p.addParameter('calculate', true);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      beam = beam@ott.beam.BscFinite(unmatched{:});
      beam.waist = p.Results.waist;
      beam.lmode = p.Results.lmode;
      beam.porder = p.Results.porder;
      beam.ellipticity = p.Results.ellipticity;
      beam.parity = p.Results.parity;
      beam.polfield = p.Results.polfield;
      beam.polbasis = p.Results.polbasis;
      beam.mapping = p.Results.mapping;
      beam.truncation_angle = p.Results.truncation_angle;
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

      beam.data = ott.bsc.LgParaxial.FromIgMode(...
          beam.waist ./ beam.wavelength, ...
          beam.lmode, beam.porder, beam.parity, beam.ellipticity, ...
          beam.polbasis, beam.polfield, ...
          'truncation_angle', beam.truncation_angle);
      beam.data = beam.data ./ beam.data.power .* beam.power;

      % Shrink Nmax to make things faster
      beam.data = beam.data.shrinkNmax('RelTol', 1e-3);
    end
  end

  methods (Hidden)
    function val = defaultVisRangeInternal(beam)
      val = [2,2] * beam.waist;
    end
  end
end

