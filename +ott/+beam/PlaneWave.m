classdef PlaneWave < ott.beam.BscInfinite
% VSWF representation of a Plane Wave beam.
% Inherits from :class:`+ott.+beam.BscInfinite`.
%
% Plane waves support smart translations, where the beam components
% are phase shifted rather than re-calculated.
%
% Properties
%   - polarisation -- (2 numeric) Polarisation in the x/y directions
%   - Nmax         -- Nmax of the stored data
%   - data         -- Internal BSC instance describing beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    polarisation    % (2 numeric) Polarisation in the x/y directions
  end

  properties (Dependent)
    Nmax            % Current Nmax of the stored data
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
      %   - initial_Nmax (numeric) -- Initial beam Nmax.  Default: ``20``.
      %     This parameter automatically grows when the beam is used but
      %     can be explicitly set for repeated use.

      p = inputParser;
      p.addOptional('polarisation', [1, 1i]);
      p.addParameter('initial_Nmax', 20);
      p.parse(varargin{:});

      beam.polarisation = p.Results.polarisation;
      beam = beam.recalculate(p.Results.initial_Nmax);
    end

    function bsc = ott.bsc.Bsc(beam, Nmax)
      % Get the BSC data for a specific Nmax
      %
      % If the current data is empty or the current Nmax is less than
      % the desired Nmax, calls :meth:`recalculate` and then returns the
      % internal data.
      %
      % Translations and rotations are applied without changing Nmax.
      %
      % Usage
      %   bsc = ott.bsc.Bsc(beam, Nmax)
      %
      % Parameters
      %   - beam (ott.beam.BscBeam) -- The beam to get data from.
      %
      %   - Nmax (numeric) -- Desired Nmax.  Default: ``0``.
      %     Resulting beam may have larger Nmax.

      if nargin < 2
        Nmax = 0;
      end

      if Nmax > beam.Nmax
        beam = beam.recalculate(Nmax);
      end

      % Apply rotations and translations (doesn't change Nmax)
      bsc = beam.data;
      bsc = bsc.rotate(beam.rotation);
      bsc = bsc.translateXyz(beam.position ./ beam.wavelength);

      % Apply polarisation
      bsc = bsc * beam.polarisation(:);
    end

    function beam = recalculate(beam, Nmax)
      % Re-calculate BSC data for specified Nmax.

      direction = repmat(beam.rotation(:, 3), 1, 2);
      poldirection = beam.rotation(:, 1:2);

      % Calculate plane waves for two components
      beam.data = ott.bsc.PlaneWave.FromDirection(...
          Nmax, direction, poldirection);

      % Clear current rotation (already applied to BSC data)
      beam.rotation = eye(3);
    end
  end

  methods % Getters/setters
    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'polarisation must be 2 element numeric vector');
      beam.polarisation = [val(1), val(2)];
    end

    function val = get.Nmax(beam)
      val = max([beam.data.Nmax]);
    end
    function beam = set.Nmax(beam, val)
      beam = beam.recalculate(val);
    end
  end
end

