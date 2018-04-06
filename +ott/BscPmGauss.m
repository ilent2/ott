classdef BscPmGauss < ott.BscPointmatch
%BscPmGauss provides HG, LG and IG beams using point matching method
%
% BscPmGauss properties:
%   type                Type of beam ('gaussian', 'lg', 'hg', or 'ig')
%   mode                Beam modes (2 or 4 element vector)
%   polarisation        Beam polarisation
%   truncation_angle    Truncation angle for beam
%   beam_offset         Offset for original beam calculation
%   waist0              Beam waist at focal plane
%
% BscPmGauss methods:
%
% This class is based on bsc_pointmatch_farfield.m and
% bsc_pointmatch_focalplane.m from ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    type               % Type of beam ('gaussian', 'lg', 'hg', or 'ig')
    mode               % Beam modes (2 or 4 element vector)
    polarisation       % Beam polarisation
    truncation_angle   % Truncation angle for beam
    beam_offset        % Offset for original beam calculation
    waist0             % Beam waist at focal plane
  end

  methods
    function beam = BscPmGauss(varargin)
      %BSCPMGAUSS construct a new IG, HG or LG gaussian beam.
      %
      % BSCPMGAUSS() constructs a new Gassian beam (LG00).
      %
      % BSCPMGAUSS(type, mode) constructs a new beam with the given type.
      % Supported types [mode]:
      %     'lg'    Laguarre-Gauss  [ radial azimuthal ]
      %     'hg'    Hermite-Gauss   [ m n ]
      %     'ig'    Ince-Gauss      [ paraxial azimuthal parity elipticity ]
      %
      % BSCPMGAUSS(..., 'Nmax') specifies the desired Nmax for the beam.
      % If omitted, Nmax is estimated using ka2nmax(k_medium*waist0).
      %
      % TODO: Documentation

      beam = beam@ott.BscPointmatch(varargin{:});
      beam.type = 'incomming';

      % Parse inputs
      p = inputParser;
      p.addOptional('type', 'lg');
      p.addOptional('mode', [ 0 0 ]);
      p.addParameter('Nmax', []);
      p.parse(varargin{:});

    end
  end
end
