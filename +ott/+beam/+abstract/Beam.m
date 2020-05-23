classdef (Abstract) Beam < ott.beam.utils.BeamInterface ...
    & matlab.mixin.Heterogeneous ...
% Base class for abstract beam descriptions (no fields).
% Inherits from :class:`ott.beam.utils.BeamInterface` and
% :class:`matlab.mixin.Hetrogeneous`.
%
% Classes that inherit from this class include a description of the
% beam but do not include methods to calculate the fields.
% These classes define casts to other beam descriptions.
%
% The class implements the :class:`BeamInterface`, which defines methods for
% calculating fields and forces.  To do this, the class casts the
% :class:`ott.beam.abstract.Beam` object to a :class:`ott.beam.Beam` object.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Beam(varargin)
      % Construct an abstract beam instance
      %
      % Usage
      %   beam = beam@ott.beam.abstract.Beam(...)
      %
      % For optional parameters, see :class:`ott.beam.properties.Beam`.

      beam = beam@ott.beam.properties.Beam(varargin{:});
    end
  end

  methods (Hidden)
    function varargout = efieldInternal(beam, varargin)
      % Cast the beam to a Beam and call method

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.efieldInternal(varargin{:});
    end

    function varargout = hfieldInternal(beam, varargin)
      % Cast the beam to a Beam and call method

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.hfieldInternal(varargin{:});
    end

    function varargout = efarfieldInternal(beam, varargin)
      % Cast the beam to a Beam and call method

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.efarfieldInternal(varargin{:});
    end

    function varargout = hfarfieldInternal(beam, varargin)
      % Cast the beam to a Beam and call method

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.hfarfieldInternal(varargin{:});
    end
  end

  methods (Static, Sealed, Access = protected)
    function default_object = getDefaultScalarElement
      % Default object for a Shape array is an empty shape
      default_object = ott.beam.abstract.Empty();
    end
  end
end

