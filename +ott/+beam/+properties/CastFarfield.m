classdef CastFarfield < ott.beam.Beam
% Base class with casts for far-field.
% Inherits from :class:`ott.beam.Beam`.
%
% Methods
%   - intensityMoment
%   - visualiseFarfield
%   - visualiseFarfieldSlice
%   - visualiseFarfieldSphere
%   - efarfieldInternal (Hidden)
%   - hfarfieldInternal (Hidden)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function varargout = visualiseFarfield(beam, varargin)
      % Create a visualisation of the beam by projecting the far-field
      % onto a plane.
      %
      % Applies a cast and calls the corresponding method of the new type.
      %
      % Usage
      %   beam.visualiseFarfield(...) displays an image of the beam
      %   in the current axes.
      %
      % See :class:`ott.beam.Beam` for parameters/usage.

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.visualiseFarfield(varargin{:});
    end

    function varargout = visualiseFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
      %
      % Applies a cast and calls the corresponding method of the new type.
      %
      % Usage
      %   beam.visualiseFarfieldSphere(...)
      %   Generate a visualisation of the far-field in the current axes.
      %
      % See :class:`ott.beam.Beam` for parameters/usage.

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.visualiseFarfieldSphere(varargin{:});
    end

    function varargout = visualiseFarfieldSlice(beam, varargin)
      % Generate a 2-D slice through the far-field
      %
      % Applies a cast and calls the corresponding method of the new type.
      %
      % Usage
      %   beam.visualiseFarfieldSlice(phi, ...)
      %   Generates a 2-D slice at angle phi around the z-axis.
      %   Plots into the current axes.
      %
      % See :class:`ott.beam.Beam` for parameters/usage.

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.visualiseFarfieldSlice(varargin{:});
    end

    function varargout = intensityMoment(beam, varargin)
      % Calculate moment of the beam intensity in the far-field
      %
      % Applies a cast and calls the corresponding method of the new type.
      %
      % Usage
      %   [moment, int] = beam.intensityMoment(...)
      %
      % See :class:`ott.beam.Beam` for parameters/usage.

      beam = ott.beam.Beam(beam);
      [varargout{1:nargout}] = beam.intensityMoment(varargin{:});
    end
  end

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
