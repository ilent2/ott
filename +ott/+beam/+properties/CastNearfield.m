classdef CastNearfield < ott.beam.Beam
% Base class with casts for near-field.
% Inherits from :class:`ott.beam.Beam`.
%
% Methods
%   - visualise
%   - efieldInternal (Hidden)
%   - hfieldInternal (Hidden)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function varargout = visualise(beam, varargin)
      % Create a visualisation of the beam
      %
      % Applies a cast and calls the corresponding method of the new type
      %
      % Usage
      %   beam.visualise(...) displays an image of the beam in the current
      %   figure window.
      %
      % See :class:`ott.beam.Beam` for parameters/usage.

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.visualise(varargin{:});
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
  end
end
