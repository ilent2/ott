classdef CastNearfield
% Base class with casts for near-field.
% Inherits from :class:`ott.beam.Beam`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

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
  end
end
