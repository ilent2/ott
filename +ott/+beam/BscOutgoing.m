classdef BscOutgoing < ott.beam.BscFinite
% Encapsulates a finite BSC beam with an outgoing basis.
% Inherits from :class:`BscFinite`.
%
% This class describes a finite BSC beam with an outgoing basis.  Fields are
% typically only valid outside the Nmax region.  The beam can be translated
% to any location outside this region.
%
% This class overrides the near-field and :class:`Bsc` cast methods to
% use outgoing VSWFs for the translations/field calculations.
%
% Additional properties/methods inherits from :class:`BscFinite`.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function beam = BscOutgoing(varargin)
      % Construct a new outgoing Finite Bsc instance
      %
      % Usage
      %   beam = BscOutgoing(data, ...)
      %
      % Optional named parameters
      %   - data (ott.bsc.Bsc) -- Initial data for beam.
      %     Default: ``ott.bsc.Bsc.empty()``.
      %
      % See base class for usage/arguments.

      beam = beam@ott.beam.BscFinite(varargin{:});
    end

    function varargout = efieldRtp(beam, varargin)
      % Calculate the electric field (in SI units).
      %
      % Uses an outgoing basis for field calculations.
      %
      % Usage
      %   [E, ...] = beam.efieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      [varargout{1:nargout}] = efieldRtp@ott.beam.BscFinite(beam, ...
          varargin{:}, 'basis', 'outgoing');
    end

    function varargout = hfieldRtp(beam, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Uses an outgoing basis for field calculations.
      %
      % Usage
      %   [H, ...] = beam.hfieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      [varargout{1:nargout}] = hfieldRtp@ott.beam.BscFinite(beam, ...
          varargin{:}, 'basis', 'outgoing');
    end
  end

  methods (Hidden)
    function bsc = translateBscInternal(beam, bsc, Nmax)
      % Applies the translation to the beam
      % Also check that the translation excludes the origin

      % Check if we have work to do
      if all(beam.position == 0)
        return;
      end

      % Check translation distance
      kanmax = ott.utils.nmax2ka(Nmax);
      kabsc = ott.utils.nmax2ka(bsc.Nmax);
      if beam.position * beam.wavenumber < kanmax + kabsc
        warning('ott:beam:BscOutgoing:small_translation', ...
          'Small translation of outgoing beam may give unexpected results');
      end

      bsc = bsc.translateXyz(beam.position ./ beam.wavelength, ...
          'Nmax', Nmax, 'basis', 'outgoing');
    end
  end
end
