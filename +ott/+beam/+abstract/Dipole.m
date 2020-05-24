classdef Dipole < ott.beam.properties.Dipole ...
    & ott.beam.abstract.Beam
% Abstract description of field from a single dipole.
% Inherits from :class:`ott.beam.properties.Dipole` and :class:`Beam`.
%
% For arrays of dipoles, see :class:`ott.beam.Dipole`.
%
% Properties
%   - polarization  -- Polarization of the dipole
%   - position      -- Position of the dipole
%   - rotation      -- Orientation of the dipole
%
% Supported casts
%   - Beam          -- Default Beam cast, uses Dipole
%   - Array         -- (Inherited from abstract.Beam)
%   - Coherent      -- Creates coherent dipole beam, Uses Dipole
%   - Incoherent    -- (Inherited from abstract.Beam)
%   - Dipole
%   - vswf.Bsc

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    polarization    % Dipole polarization
  end

  methods
    function beam = Dipole(varargin)
      % Construct a new dipole representation
      %
      % Usage
      %   beam = Dipole(polarization, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - polarization (3 numeric) -- Polarization of dipole [x;y;z].

      beam = beam@ott.beam.properties.Dipole(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Construct a new Dipole instance
      beam = ott.beam.Dipole(varargin{:});
    end

    function beam = ott.beam.Dipole(beam, varargin)
      % Construct a new Dipole instance
      %
      % Arrays of dipoles are assumed to be coherent.  Use the
      % Array or Incoherent classes if this is not the case.

      assert(isa(beam, 'ott.beam.abstract.Dipole'), ...
          'First argument must be abstract.Dipole');

      % Handle coherent beam arrays
      position = [beam.position];
      polarization = reshape([beam.polarization], [], 1);

      % TODO: Other properties
      beam = ott.beam.Dipole(position, polarization);
    end

    function beam = ott.beam.vswf.Bsc(beam, varargin)

      assert(isa(beam, 'ott.beam.abstract.Dipole'), ...
          'First argument must be abstract.Dipole');

      % TODO: Support for arrays of dipoles
      % TODO: Thing about this, do it smartly...
      beam = ott.beam.vswf.Bsc(a, b, 'basis', 'outgoing');
    end

    function beam = ott.beam.Coherent(varargin)
      % Convert to Dipole
      beam = ott.beam.Dipole(varargin{:});
    end
  end

  methods % Getters/setters
    function beam = set.polarization(beam, val)
      assert(isnumeric(val) && isvector(val) && numel(val) == 3, ...
          'polarization must be 3-element numeric vector');
      beam.polarization = val;
    end
  end
end
