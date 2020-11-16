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
      % Optional named parameters
      %   - truncation_angle (numeric) -- Angle to truncate the beam.
      %     Defaults to match the beam NA: ``asin(NA/index_medium)``.
      %
      % For other parameters, see class constructor.

      p = inputParser;
      p.addOptional('mmode', 0, @isnumeric);
      p.addOptional('nmode', 0, @isnumeric);
      p.addParameter('mapping', 'sin');
      p.addParameter('index_medium', 1.0);
      p.addParameter('truncation_angle', []);
      p.addParameter('calculate', true);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      paraxial_order = p.Results.mmode + p.Results.nmode;

      % Set the default truncation angle to the specified NA
      truncation_angle = p.Results.truncation_angle;
      if isempty(truncation_angle)
        truncation_angle = ott.utils.na2angle(Na, p.Results.index_medium);
      end

      % Construct beam (parsing remaining parameters)
      beam = ott.beam.HermiteGaussian(0, ...
          'truncation_angle', truncation_angle, ...
          varargin{:}, 'calculate', false);
        
      % Calculate and set waist
      waist = ott.bsc.LgParaxial.WaistFromNa(Na, p.Results.index_medium, ...
          paraxial_order, p.Results.mapping);
      beam.waist = waist * beam.wavelength;
      
      % Calculate beam data if requested
      if p.Results.calculate
        beam = beam.recalculate([]);
      end
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
      %   - truncation_angle (numeric) -- Adds a hard edge truncation
      %     to the beam in the far-field.  Default: ``pi/2`` (i.e., no
      %     truncation added).
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
      p.addParameter('truncation_angle', pi/2);
      p.addParameter('calculate', true);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      beam = beam@ott.beam.BscFinite(unmatched{:});
      beam.waist = p.Results.waist;
      beam.mmode = p.Results.mmode;
      beam.nmode = p.Results.nmode;
      beam.polfield = p.Results.polfield;
      beam.polbasis = p.Results.polbasis;
      beam.mapping = p.Results.mapping;
      beam.truncation_angle = p.Results.truncation_angle;

      if p.Results.calculate
        beam = beam.recalculate([]);
      end
    end

    function varargout = recalculate(beam, ~, varargin)
      % Re-calculate data for beam.
      %
      % This function can be called to pre-compute the beam data for
      % the current beam parameters.  It is automatically called when
      % the beam is used, however, calling the function explicitly can
      % speed up run-time.
      %
      % Usage
      %   [beam, ...] = beam.recalcualte(Nmax, ...)
      %
      % Parameters
      %   - Nmax (numeric) -- Parameter is ignored.
      %
      % Unmatched arguments passed to, and additional results returned
      % from :meth:`LgParaxial.FromHgMode`.

      ott.utils.nargoutCheck(beam, nargout);

      [beam.data, weights, varargout{2:nargout}] = ...
          ott.bsc.LgParaxial.FromHgMode(...
          beam.waist ./ beam.wavelength, ...
          beam.mmode, beam.nmode, beam.polbasis, beam.polfield, ...
          'truncation_angle', beam.truncation_angle, varargin{:});
      beam.data = beam.data * weights;

      % Shrink Nmax to make things faster
      beam.data = beam.data.shrinkNmax('RelTol', beam.getSetShrinkNmaxRelTol);

      % Normalise beam power
      beam.data.power = 1;

      varargout{1} = beam;
    end
  end

  methods (Hidden)
    function val = defaultVisRangeInternal(beam)
      val = [2,2] * beam.waist;
    end
  end
end

