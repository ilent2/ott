classdef BscInfinite < ott.beam.BscBeam
% Describes a beam with infinite spatial extent.
% Inherits from :class:`BscBeam`.
%
% This class is useful for describing plane waves, annular beams,
% and other beams with an infinite spatial extent.  Values are only
% valid within the Nmax region, requesting values outside this region
% requires the beam to be re-calculated (or, for plane waves/annular
% beams, translated).
%
% This class overloads the field calculation functions and requests
% a larger Nmax whenever points outside the valid range are requested.
% Far-field functions raise a warning that fields may not look as expected.
%
% For methods/properties, see :class:`BscBeam`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = BscInfinite(varargin)
      % Construct a new BSC infinite beam.
      %
      % Usage
      %   beam = BscInfinite(data, ...)
      %
      % Optional named parameters
      %   - data (ott.bsc.Bsc) -- Initial data for beam or empty.
      %     Default: ``ott.bsc.Bsc.empty()``.
      %
      %   - index_medium (numeric) -- Refractive index of the medium.
      %     Default: ``1.0``.
      %
      %   - omega (numeric) -- Optical angular frequency [Hz].
      %     Default: ``3e8/1064e-9*2*pi`` (i.e., default vacuum
      %     wavelength is 1064 nm).

      beam = beam@ott.beam.BscBeam(varargin{:});
    end

    %
    % Updated field calculation functions
    %

    % efield, hfield, and ehfield are all fine (they use Rtp methods)
    % hfarfield and ehfarfield both call efarfield

    function varargout = efarfield(beam, varargin)
      % Calculate the electric field (in SI units)
      %
      % Raises a warning that far-fields may be inaccurate when using
      % a finite VSWF representation of an infinite beam.
      %
      % Usage
      %   [E, H, ...] = beam.efarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      warning('ott:beam:BscInfinite:farfield_is_finite', ...
          'Far-fields may be inaccurate for infinite beams');

      [varargout{1:nargout}] = efarfield@ott.beam.BscBeam(beam, varargin{:});
    end

    function varargout = ehfieldRtp(beam, rtp, varargin)
      % Calculate E and H fields in spherical coordinates. (SI units)
      %
      % If points are outside the Nmax range, re-calculate beam data.
      %
      % Usage
      %   [E, H, data] = beam.ehfieldRtp(rtp, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r; theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      % Calculate required Nmax
      oNmax = ott.utils.ka2nmax(max(rtp(1, :)).*beam.wavenumber);

      % Create new beam with updated data and calculate fields
      beam = ott.beam.BscBeam(ott.bsc.Bsc(beam, oNmax));
      [varargout{1:nargout}] = beam.ehfieldRtp(rtp, varargin{:});
    end

    function varargout = efieldRtp(beam, rtp, varargin)
      % Calculate the electric field (in SI units)
      %
      % If points are outside the Nmax range, re-calculate beam data.
      %
      % Usage
      %   [E, ...] = beam.efieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      % Calculate required Nmax
      oNmax = ott.utils.ka2nmax(max(rtp(1, :)).*beam.wavenumber);

      % Create new beam with updated data and calculate fields
      beam = ott.beam.BscBeam(ott.bsc.Bsc(beam, oNmax));
      [varargout{1:nargout}] = beam.efieldRtp(rtp, varargin{:});
    end

    function varargout = hfieldRtp(beam, rtp, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % If points are outside the Nmax range, re-calculate beam data.
      %
      % Usage
      %   [H, ...] = beam.hfieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      % Calculate required Nmax
      oNmax = ott.utils.ka2nmax(max(rtp(1, :)).*beam.wavenumber);

      % Create new beam with updated data and calculate fields
      beam = ott.beam.BscBeam(ott.bsc.Bsc(beam, oNmax));
      [varargout{1:nargout}] = beam.hfieldRtp(rtp, varargin{:});
    end
  end
end
