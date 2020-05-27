classdef FocussedTopHat < ott.beam.properties.FocussedTopHat ...
    & ott.beam.properties.MaterialBeam ...
    & ott.beam.abstract.Beam
% Abstract description of a focussed top-hat beam.
% Inherits from :class:`ott.beam.properties.FocussedTopHat` and :class:`Beam`.
%
% Properties
%   - angle
%
% TODO: The Beam cast looks wrong, we should have a full
%     beam not an abstract
%
% Supported casts
%   - Beam                -- Uses abstract.MaskedFarfield
%   - Ray
%   - abstract.InterpFarfield   -- (Inherited)
%   - abstract.MaskedFarfield
%   - abstract.MaskedNearfield3d
%   - vswf.Bsc            -- (Inherited) Uses vswf.FarfieldPm
%   - vswf.Pointmatch     -- (Inherited) Uses vswf.FarfieldPm
%   - vswf.FarfieldPm     -- (Inherited) Uses Beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.properties.FocussedTopHat.likeProperties(other, args);
      args = ott.beam.abstract.Beam.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = UniformFarfield.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.FocussedTopHat.likeProperties(other, varargin);
      beam = ott.beam.abstract.FocussedTopHat(args{:});
    end
  end

  methods
    function beam = FocussedTopHat(varargin)
      % Constructs a new representation of a focussed top hat beam
      %
      % Usage
      %   beam = FocussedTopHat(angle, ...)

      beam = beam@ott.beam.properties.FocussedTopHat(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Cast the beam to a MaskedFarfield
      beam = ott.beam.MaskedFarfield(varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Cast to collection of rays

      % TODO: Multiple beams and type check

      theta = linspace(0, beam.angle, 50);
      phi = linspace(0, 2*pi, 100);
      radius = 100;
      [R, T, P] = meshgrid(radius, theta, phi);
      vrtp = [0*R(:), ones(size(T(:))), 0*P(:)].';
      [vxyz, xyz] = ott.utils.rtp2xyz(vrtp, [R(:), T(:), P(:)].');

      beam = ott.beam.Ray('origin', xyz, 'dierction', -xyz, ...
          'polarisation1', vxyz, varargin{:});
    end

    function beam = ott.beam.abstract.MaskedFarfield(beam, varargin)
      % Construct a masked far-field beam
      %
      % Usage
      %   beam = MaskedFarfield(beam, masked_beam, ...)
      %   For additional arguments see :class:`MaskedFarfield`.

      % TODO: Multiple beams and type check

      masked_beam = ott.beam.abstract.UniformFarfield();

      beam = ott.beam.abstract.MaskedFarfield.TopHat(beam.angle, ...
          masked_beam, varargin{:});
    end

    function beam = ott.beam.abstract.MaskedNearfield3d(beam, varargin)
      % Construct a masked near-field beam
      %
      % Usage
      %   beam = MaskedNearfield3d(beam, masked_beam, ...)
      %   For additional arguments see :class:`MaskedNearfield3d`.

      % TODO: Multiple beams and type check

      rtp_mask = @(rtp) rtp(2, :) < beam.angle;
      xyz_mask = @(xyz) rtp_mask(ott.utils.xyz2rtp(xyz));

      masked_beam = ott.beam.abstract.UniformFarfield();

      beam = ott.beam.abstract.MaskedNearfield3d(xyz_mask, ...
          masked_beam, varargin{:});
    end
  end
end
