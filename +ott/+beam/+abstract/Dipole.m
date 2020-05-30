classdef Dipole < ott.beam.properties.Dipole ...
    & ott.beam.abstract.CastBoth
% Abstract description of field from a single dipole.
% Inherits from :class:`ott.beam.properties.Dipole` and :class:`CastBoth`.
%
% For arrays of dipoles, see :class:`ott.beam.Dipole`.
%
% Static methods
%   - like            -- Create a beam like another
%
% Properties
%   - location      -- Alias for dipole position
%   - polarization  -- Polarization of the dipole
%   - position      -- Position of the dipole
%   - rotation      -- Orientation of the dipole
%
% Casts
%   - Beam          -- Uses Dipole
%   - Coherent      -- Uses Dipole (no need for beam.Array)
%   - Dipole
%   - vswf.Bsc      -- Uses vswf.Dipole
%   - vswf.Dipole
%
% See base for other supported casts.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    polarization    % Dipole polarization
  end

  properties (Dependent, SetAccess=protected)
    location        % Alias for position
  end

  properties
    power           % Power of dipole
  end

  methods (Static)
    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Dipole.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.Dipole.likeProperties(other, varargin);
      beam = ott.beam.abstract.Dipole(args{:});
    end
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
      %
      %   - location (3 numeric) -- Location of dipole (alias for location).
      %   - position (3 numeric) -- Position of dipole.  Default: ``[0,0,0]``.

      p = inputParser;
      p.addOptional('polarization', [], @isnumeric);
      p.addParameter('position', [0;0;0]);
      p.addParameter('location', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if ~any(strcmpi(p.UsingDefaults, 'location')) ...
          && ~any(strcmpi(p.UsingDefaults, 'position'))
        warning('ott:beam:abstract:Dipole:position_and_location', ...
            'Both position and location specified, ignoring position');
        location = p.Results.location;
      elseif ~any(strcmpi(p.UsingDefaults, 'location'))
        location = p.Results.location;
      else
        location = p.Results.position;
      end

      beam = beam@ott.beam.properties.Dipole(unmatched{:}, ...
          'location', location, 'polarization', p.Results.polarization);
    end

    function beam = ott.beam.Beam(varargin)
      % Construct a new Dipole instance
      beam = ott.beam.Dipole(varargin{:});
    end

    function beam = ott.beam.Coherent(varargin)
      % Convert to Dipole
      beam = ott.beam.Dipole(varargin{:});
    end

    function beam = ott.beam.Dipole(beam, varargin)
      % Construct a new Dipole instance
      beam = castArrayHelper(@ott.beam.vswf.Dipole.like, beam, varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to vswf.Dipole
      beam = ott.beam.vswf.Dipole(varargin{:});
    end

    function beam = ott.beam.vswf.Dipole(beam, varargin)
      % Cast to vswf.Dipole
      beam = castArrayHelper(@ott.beam.vswf.Dipole.like, beam, varargin{:});
    end
  end

  methods (Access=protected)
    function beam = castArrayHelper(cast, beam, varargin)
      % Construct a new Dipole instance
      %
      % Arrays of dipoles are assumed to be coherent.  Use the
      % Array or Incoherent classes if this is not the case.

      assert(isa(beam, 'ott.beam.abstract.Dipole'), ...
          'First argument must be abstract.Dipole');

      ott.utils.nargoutCheck(beam, nargout);

      % Handle coherent beam arrays
      arg_location = [beam.position];
      arg_polarization = reshape([beam.polarization], [], 1);

      beam = cast(beam, 'position', [0;0;0], 'rotation', eye(3), ...
          'location', arg_location, 'polarization', arg_polarization, ...
          varargin{:});
    end
  end

  methods % Getters/setters
    function beam = set.polarization(beam, val)
      assert(isnumeric(val) && isvector(val) && numel(val) == 3, ...
          'polarization must be 3-element numeric vector');
      beam.polarization = val;
    end

    function l = get.location(beam)
      l = beam.position;
    end
    function beam = set.location(beam, val)
      beam.position = val;
    end

    function p = get.power(~)
      error('Not yet implemented');
    end
  end
end
