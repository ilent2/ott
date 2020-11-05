classdef BscBeam < ott.beam.ArrayType
% Beam class encapsulating a BSC instance.
% Inherits from :class:`+ott.+beam.ArrayType`.
%
% This class stores an internal :class:`+ott.+bsc.Bsc` instance which it
% uses for calculating fields.
%
% Methods in this class assume the beam is a regular beam (i.e., not
% the outgoing fields from a scattered beam).
%
% Properties
%   - data          -- Internal BSC instance describing beam
%   - defaultVisRange -- Default range for near-field visualisations
%
% Methods
%   - recalculate   -- Can be overloaded by sub-classes to update data
%   - efield, hfield, ... -- Calculate fields in SI units, uses the
%     :class:`ott.bsc.Bsc` methods internally.
%   - force, torque, spin -- Calculate force/torque/spin in SI units
%
% Supported casts
%   - ott.bsc.Bsc   -- Get the BSC data after applying transformations
%
% Additional properties/methods inherited from :class:`Beam`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    data      % Internal BSC instance describing beam
  end

  properties (Dependent)
    defaultVisRange     % Chosen based on beam NA
  end

  methods
    function beam = BscBeam(varargin)
      % Construct a Bsc beam instance.
      %
      % Usage
      %   beam = BscBeam(data, ...)
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

      p = inputParser;
      p.addOptional('data', ott.bsc.Bsc.empty(), @(x) isa(x, 'ott.bsc.Bsc'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.ArrayType(unmatched{:});
      beam.data = p.Results.data;
    end

    function bsc = ott.bsc.Bsc(beam, Nmax)
      % Get the BSC data for a specific Nmax.
      %
      % If the current data is empty or the current Nmax is less than
      % the desired Nmax, calls :meth:`recalculate` and then returns the
      % internal data.
      %
      % If position and rotation have been set, rotations and translations
      % are applied to the beam (rotations first, translations second).
      % Translations of this beam type require a larger internal Nmax,
      % if the Nmax is not large enough for the translation,
      % :meth:`recalculate` is called.
      %
      % Usage
      %   bsc = ott.bsc.Bsc(beam, Nmax)
      %
      % Parameters
      %   - beam (ott.beam.BscBeam) -- The beam to get data from.
      %
      %   - Nmax (numeric) -- Desired Nmax.  The default is either 0
      %     or the current beam Nmax minus the Nmax required for any
      %     translations, whichever is greater.

      % Calculate Nmax required after translation
      % This assumes the beam is a regular beam
      % Also assumes generic translation
      dist = vecnorm(beam.position);
      tNmax = ott.utils.ka2nmax(beam.wavenumber*dist);

      if nargin < 2
        Nmax = max(0, max([0, beam.data.Nmax]) - tNmax);
      end
      assert(isnumeric(Nmax) && isscalar(Nmax), ...
        'Nmax must be numeric scalar');

      if isempty(beam.data) || max([beam.data.Nmax]) < (Nmax + tNmax)
        beam = beam.recalculate(Nmax + tNmax);
      end

      bsc = beam.data;

      % Rotate the beam
      bsc = bsc.rotate(beam.rotation);

      % Translate the beam
      % Assumes beam is a regular beam
      bsc = bsc.translateXyz(beam.position ./ beam.wavelength, 'Nmax', Nmax);

      % Combine beams now if coherent
      if strcmpi(beam.arrayType, 'coherent')
        bsc = sum(bsc);
      end
    end

    function beam = recalculate(beam, Nmax)
      % Re-calculate data for a specified Nmax
      %
      % This implementation simply raises an error when called.
      % This function should be overloaded by sub-classes of BscBeam
      % to implement the desired behaviour.
      %
      % Usage
      %   beam = beam.recalcualte(Nmax)
      %
      % Parameters
      %   - Nmax (numeric) -- Desired Nmax for beam.

      error('ott:beam:BscBeam:recalculate_not_implemented', ...
          'recalculate is not implemented for these type of beam');
    end

    %
    % Field calculation functions
    %

    function varargout = efarfield(beam, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, H, ...] = beam.efarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      bsc = ott.bsc.Bsc(beam);
      [varargout{1:nargout}] = bsc.efarfield(varargin{:});

      if nargout >= 1
        varargout{1} = varargout{1} ./ beam.wavenumber;
      end
    end

    function [H, vswfData] = hfarfield(beam, rtp, varargin)
      % Calculate H far-field (uses :meth:`efarfield`)
      %
      % Usage
      %   [H, data] = beam.hfarfield(rtp, ...)
      %
      % See :meth:`efarfield` for further details.

      [~, H, vswfData] = beam.ehfarfield(rtp, varargin{:});
    end

    function [E, H, data] = ehfarfield(beam, rtp, varargin)
      % Calculate E and H far-fields
      %
      % Usage
      %   [E, H, data] = beam.ehfarfield(rtp, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      % Calculate E as normal
      [E, data] = beam.efarfield(rtp, varargin{:});

      % Swap theta and phi components
      H = ott.utils.FieldVector(rtp, ...
          -1i .* E.vrtp([1, 3, 2], :), 'spherical');
      H = H ./ beam.impedance;
    end

    function varargout = efield(beam, xyz, varargin)
      % Calculate E field from Cartesian coordinates.
      %
      % Usage
      %   [E, data] = beam.efield(xyz, ...)
      %   E is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.efieldRtp(rtp, varargin{:});
    end

    function varargout = hfield(beam, xyz, varargin)
      % Calculate H field from Cartesian coordinates.
      %
      % Usage
      %   [H, data] = beam.hfield(xyz, ...)
      %   H is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.hfieldRtp(rtp, varargin{:});
    end

    function varargout = ehfield(beam, xyz, varargin)
      % Calculate E and H fields from Cartesian coordinates.
      %
      % Usage
      %   [E, H, data] = beam.ehfield(xyz, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      % This overloads the default in order to avoid two xyz2rtp casts

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.ehfieldRtp(rtp, varargin{:});
    end

    function varargout = efieldRtp(beam, rtp, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, ...] = beam.efieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      bsc = ott.bsc.Bsc(beam);
      [varargout{1:nargout}] = bsc.efieldRtp(...
          rtp./[beam.wavelength;1;1], varargin{:});

      if nargout >= 1
        varargout{1} = varargout{1} ./ beam.wavenumber;
      end
    end

    function varargout = hfieldRtp(beam, rtp, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [H, ...] = beam.hfieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      bsc = ott.bsc.Bsc(beam);
      [varargout{1:nargout}] = bsc.hfieldRtp(...
          rtp./[beam.wavelength;1;1], varargin{:});

      if nargout >= 1
        varargout{1} = varargout{1} ./ beam.wavenumber ./ beam.impedance;
      end
    end

    function varargout = force(ibeam, sbeam, varargin)
      % Calculate force (in Newtons)
      %
      % Usage
      %   force = incident_beam.force(scattered_beam) -- Calculates
      %   difference in momentum between two beams.
      %
      %   force = incident_beam.force(particle, ...) -- First calculates
      %   scattering, then calculates force from difference in momentum.
      %   Additional parameters are passed to :meth:`scatter`.

      % Calculate scattering if required
      if ~isa(sbeam, 'ott.beam.Beam')
        [sbeam, ibeam] = ibeam.scatter(sbeam, varargin{:});
      end

      % Get data from beams
      ibsc = ott.bsc.Bsc(ibeam);
      sbsc = ott.bsc.Bsc(sbeam);

      % Calculate force using internal methods
      [varargout{1:nargout}] = ibsc.force(sbsc);

      % Convert units of output to SI
      for ii = 1:nargout
        varargout{ii} = varargout{ii} ./ ibeam.speed;
      end
    end

    function varargout = torque(ibeam, sbeam)
      % Calculate torque (in Newton meters)
      %
      % Usage
      %   torque = incident_beam.torque(scattered_beam) -- Calculates
      %   difference in momentum between two beams.
      %
      %   torque = incident_beam.torque(particle, ...) -- First calculates
      %   scattering, then calculates torque from difference in momentum.
      %   Additional parameters are passed to :meth:`scatter`.

      % Calculate scattering if required
      if ~isa(sbeam, 'ott.beam.Beam')
        [sbeam, ibeam] = ibeam.scatter(sbeam, varargin{:});
      end

      % Get data from beams
      ibsc = ott.bsc.Bsc(ibeam);
      sbsc = ott.bsc.Bsc(sbeam);

      % Calculate torque using internal methods
      [varargout{1:nargout}] = ibsc.torque(sbsc);

      % Convert units of output to SI
      for ii = 1:nargout
        varargout{ii} = varargout{ii} ./ ibeam.omega;
      end
    end

    function varargout = spin(ibeam, sbeam)
      % Calculate spin (in Newton meters)
      %
      % Usage
      %   spin = incident_beam.spin(scattered_beam) -- Calculates
      %   difference in momentum between two beams.
      %
      %   spin = incident_beam.spin(particle, ...) -- First calculates
      %   scattering, then calculates spin from difference in momentum.
      %   Additional parameters are passed to :meth:`scatter`.

      % Calculate scattering if required
      if ~isa(sbeam, 'ott.beam.Beam')
        [sbeam, ibeam] = ibeam.scatter(sbeam, varargin{:});
      end

      % Get data from beams
      ibsc = ott.bsc.Bsc(ibeam);
      sbsc = ott.bsc.Bsc(sbeam);

      % Calculate spin using internal methods
      [varargout{1:nargout}] = ibsc.spin(sbsc);

      % Convert units of output to SI
      for ii = 1:nargout
        varargout{ii} = varargout{ii} ./ ibeam.omega;
      end
    end
  end

  methods % Getters/setters
    function beam = set.data(beam, val)
      if isempty(val)
        val = ott.bsc.Bsc.empty();
      end
      assert(isa(val, 'ott.bsc.Bsc'), ...
          'data must be a ott.bsc.Bsc instance or empty');
      beam.data = val;
    end

    function val = get.defaultVisRange(beam)
      val = [1,1] .* ott.utils.nmax2ka(...
          max([1, beam.data.Nmax]))./beam.wavenumber ./ sqrt(2);
    end
  end
end
