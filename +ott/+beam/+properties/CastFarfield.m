classdef CastFarfield
% Base class with casts for far-field.
% Inherits from :class:`ott.beam.Beam`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Hidden)
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
end
