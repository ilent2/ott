classdef FocussedTopHat < ott.beam.properties.FocussedTopHat ...
    & ott.beam.abstract.Beam
% Abstract description of a focussed top-hat beam.
% Inherits from :class:`ott.beam.properties.FocussedTopHat` and :class:`Beam`.
%
% Properties
%   - angle
%
% Supported casts
%   - Beam                -- Uses vswf.FarfieldPm
%   - Ray
%   - abstract.InterpFarfield
%   - abstract.MaskedFarfield
%   - abstract.MaskedNearfield3d
%   - vswf.Bsc            -- Uses vswf.FarfieldPm
%   - vswf.Pointmatch     -- Uses vswf.FarfieldPm
%   - vswf.FarfieldPm

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = FocussedTopHat(varargin)
      % Constructs a new representation of a focussed top hat beam
      %
      % Usage
      %   beam = FocussedTopHat(angle, ...)

      beam = beam@ott.beam.properties.FocussedTopHat(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to FarfieldPm
      beam = ott.beam.vswf.FarfieldPm(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to FarfieldPm
      beam = ott.beam.vswf.FarfieldPm(varargin{:});
    end

    function beam = ott.beam.vswf.Pointmatch(varargin)
      % Cast to FarfieldPm
      beam = ott.beam.vswf.FarfieldPm(varargin{:});
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

    function newbeam = ott.beam.abstract.InterpFarfield(beam, varargin)
      % Construct a far-field interpolated beam

      % TODO: Multiple beams and type check

      newbeam = ott.beam.abstract.MaskedFarfield(beam);

      % Calculate far-field pattern
      % TODO: Might be smarter choice
      theta = linspace(0, pi, 100);
      phi = linspace(0, 2*pi, 100);
      [T, P] = meshgrid(theta, phi);
      tp = [T(:), P(:)].';
      [Etp, Htp] = newbeam.ehfarfield(tp);

      newbeam = ott.beam.abstract.InterpFarfield(tp, Etp, Htp, varargin{:});
    end

    function beam = ott.beam.vswf.FarfieldPm(beam, varargin)
      % First cast to InterpFarfield then to FarfieldPm

      % TODO: Multiple beams and type check

      beam = ott.beam.abstract.InterpFarfield(beam);
      beam = ott.beam.vswf.FarfieldPm(beam, varargin{:});
    end
  end
end
