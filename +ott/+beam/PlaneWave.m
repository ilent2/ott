classdef PlaneWave < ott.beam.BscInfinite
% VSWF representation of a Plane Wave beam.
% Inherits from :class:`+ott.+beam.BscInfinite`.
%
% Properties
%   - polarisation -- (2 numeric) Polarisation in the x/y directions
%   - data      -- Internal BSC instance describing beam
%
% Methods
%   - getData         -- Get data for specific Nmax
%   - recalculate     -- Update the internal data for new Nmax

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    polarisation    % (2 numeric) Polarisation in the x/y directions
  end

  properties (Hidden)
    polarisationInternal
  end

  methods
    function beam = PlaneWave(varargin)
      % Construct a new plane wave beam VSWF representation.
      %
      % Usage
      %   beam = ott.beam.PlaneWave(polarisation, ...)
      %
      % Optional named arguments
      %   - polarisation (2 numeric) -- Polarisation in the x/y direction.
      %     Polarisation directions correspond to the first two columns
      %     of the ``rotation`` matrix.  When the beam is not rotated,
      %     these are the x/y directions (for a beam propagating in z).
      %     Default: ``[1, 1i]``.
      %
      %   - initial_Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('polarisation', [1, 1i]);
      p.addParameter('initial_Nmax', 0);
      p.parse(varargin{:});

      beam.polarisation = p.Results.polarisation;
      beam = beam.recalculate(p.Results.initial_Nmax);
    end
  end

  methods (Hidden)
    function [data, vswfData] = recalculateInternal(beam, Nmax, vswfData)
      % Re-calculate BSC data for specified Nmax.

      direction = repmat(beam.rotation(:, 3), 1, 2);
      poldirection = beam.rotation(:, 1:2);

      % Calculate plane waves for two components
      [data, vswfData] = ott.bsc.PlaneWave.FromDirection(...
          Nmax, direction, poldirection, 'data', vswfData);

      % Apply polarisation
      data = data * beam.polarisation(:);
    end
  end

  methods % Getters/setters
    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'polarisation must be 2 element numeric vector');
      beam.polarisationInternal = [val(1), val(2)];
      beam.data = [];
    end
    function val = get.polarisation(beam)
      val = beam.polarisationInternal;
    end
  end
end

