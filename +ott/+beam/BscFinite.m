classdef BscFinite < ott.beam.BscBeam
% A beam represented by a finite VSWF expansion.
% Inherits from :class:`BscBeam`.
%
% This class describes beams which can be represented using a finite VSWF
% expansion. The class stores a :class:`Bsc` instance internally.
% BSC coefficients at any other location can be found by applying a
% translation to the beam data.
%
% Far-fields are calculated without applying a translation to the BSC
% data, instead the fields are calculated at the origin and phase shifted.
%
% As with the :class:`BscBeam` class, this class assumes a regular beam.
%
% Properties
%   - power         -- Power applied to the beam in ``ott.bsc.bsc``.
%
% Supported casts
%   - ott.bsc.Bsc   -- Get the BSC data after applying transformations
%
% Additional properties/methods inherited from :class:`BscBeam`.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    power       % Beam power [W]
  end

  methods (Static)
    function val = getSetShrinkNmaxRelTol(val)
      % Get or set the shringNmax RelTol value for BscFinite beams.
      %
      % This parameter is used by the re-calculate methods to shrink
      % the generated beam and speed up simulations.  The default
      % parameter is ``2e-2`` but it can be changed with
      %
      %   ott.beam.BscFinite.getSetShrinkNmaxRelTol(val)

      persistent RelTol
      if isempty(RelTol)
        RelTol = 2e-2;
      end

      if nargin == 1
        RelTol = val;
      else
        val = RelTol;
      end
    end
  end

  methods
    function beam = BscFinite(varargin)
      % Construct a new Finite Bsc instance
      %
      % Usage
      %   beam = BscFinite(data, ...)
      %
      % Optional named parameters
      %   - data (ott.bsc.Bsc) -- Initial data for beam.
      %     Default: ``ott.bsc.Bsc.empty()``.
      %
      %   - index_medium (numeric) -- Refractive index of the medium.
      %     Default: ``1.0``.
      %
      %   - omega (numeric) -- Optical angular frequency [Hz].
      %     Default: ``3e8/1064e-9*2*pi`` (i.e., default vacuum
      %     wavelength is 1064 nm).
      %
      %   - power (numeric) -- Beam power [W].  Applied to the beam
      %     whenever ``ott.bsc.Bsc`` is called (i.e., when the beam is used).
      %     Default: ``1``.

      p = inputParser;
      p.addOptional('data', ott.bsc.Bsc.empty(), @(x) isa(x, 'ott.bsc.Bsc'));
      p.addParameter('power', 1, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.BscBeam(p.Results.data, unmatched{:});
      beam.power = p.Results.power;
    end

    function bsc = ott.bsc.Bsc(beam, Nmax)
      % Get the beam data for the specified Nmax.
      %
      % If the current data is empty, calls :meth:`recalculate`.
      %
      % If the position and rotation have been set, rotations and translations
      % are applied to the beam (rotations first, translations second).
      % The translated beam will have the specified Nmax.
      %
      % Usage
      %   bsc = ott.bsc.Bsc(beam, Nmax)
      %
      % Parameters
      %   - beam (ott.beam.BscFinite) -- The beam to get data from.
      %
      %   - Nmax (numeric) -- Desired Nmax.  This parameter is ignored
      %     except when a translation is applied, in which case the
      %     translated beam has this Nmax.  Default: current Nmax or 0.

      if nargin < 2
        Nmax = max([0, beam.data.Nmax]);
      end

      if isempty(beam.data)
        beam = beam.recalculate(Nmax);
      end

      bsc = beam.data;

      % Rotate the beam
      bsc = bsc.rotate(beam.rotation);

      % Translate the beam
      bsc = beam.translateBscInternal(bsc, Nmax);

      % Apply power to beam
      bsc.power = bsc.power * beam.power;

      % Apply scale
      bsc = bsc * beam.scale;
    end
    
    %
    % Field calculation functions
    %

    function varargout = efarfield(beam, varargin)
      % Calculate the electric field (in SI units)
      %
      % Calculates the fields in the local reference frame and then
      % applies phase/rotations to give the global far-field.
      %
      % Usage
      %   [E, H, ...] = beam.efarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.
      
      % Get rotation/position and clear beam properties
      % TODO: Rotations
      scat_position = beam.position;
      beam.position = [0;0;0];
      
      % Defer to base for field calculation
      [varargout{1:nargout}] = efarfield@ott.beam.BscBeam(beam, varargin{:});

      if nargout >= 1
        % Apply phase shifts to far-fields
        varargout{1} = beam.translateFarfields(varargout{1}, -scat_position);
      end
    end
  end

  methods (Hidden)
    function bsc = translateBscInternal(beam, bsc, Nmax)
      % Applies the translation to the beam shape coefficients
      % Can be overloaded by sub-classes to change default behaviour

      bsc = bsc.translateXyz(-beam.position ./ beam.wavelength, ...
          'Nmax', Nmax, 'basis', 'regular');
    end
  end

  methods % Getters/setters
    function beam = set.power(beam, val)
      assert(isnumeric(val) && isscalar(val) && val >= 0, ...
          'power must be positive numeric scalar');
      beam.power = val;
    end
  end
end
