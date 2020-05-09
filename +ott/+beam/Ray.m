classdef Ray < ott.beam.PlaneWave
% Specialisation of a PlaneWave for Ray-optics beams
% Inherits from :class:`PlaneWave`.
%
% The main differences between :class:`PlaneWave` and :class:`Ray` are
% finite power and more ray-orientated defaults for methods.
%
% Properties
%   - origin        -- Ray origins, 3xN array (default [0;0;0])
%   - direction     -- Direction of propagation (3xN Cartesian)
%   - field         -- Field parallel and perpendicular to polarisation
%   - polarisation  -- Primary polarisation direction
%
% Dependent properties
%   - power         -- Power of the beam, calculated from the field property

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function beam = empty(varargin)
      % Construct an emtpy beam array
      %
      % Usage
      %   beam = ott.beam.Ray.empty(...)
      %
      % Additional parameters are passed to the constructor.
      
      empt = zeros(3, 0);
      beam = ott.beam.Ray('direction', empt, 'polarisation', empt, ...
        'field', empt(1, :), 'origin', empt, varargin{:});
    end
  end

  methods
    function beam = Ray(varargin)
      % Construct a new geometric optics ray instance
      %
      % Usage
      %   beam = Ray(...)
      %
      % Named arguments
      %   - direction (3xN numeric) -- direction vectors (Cartesian)
      %     Default: ``[0;0;1]``.
      %
      %   - polarisation (3xN numeric) -- polarisation vectors (Cartesian)
      %     Default: ``[1;0;0]``.
      %
      %   - field (1xN|2xN numeric) -- Field vectors parallel and
      %     (optionally) perpendicular to the polarisation direction.
      %     Allows for 0 intensity with finite polarisation direction.
      %     Default: ``1``.
      %
      %   - origin (3xN numeric) -- Origin of plane waves.
      %     Default: ``[0;0;0]``.
      %
      %   - vector (ott.utils.Vector) -- Vector describing origin and
      %     direction of the Ray.  Incompatible with `direction` and
      %     `origin`.  Default: ``[]``.

      % Call the base constructor (must be outside if statement)
      beam = beam@ott.beam.PlaneWave(varargin{:});
    end

    function ray = focus(ray, location)
      % Focus rays to a point
      %
      % Rotates each ray to face the specified location.
      % For this method to do anything sensible, the origins of
      % the origin rays must be set.
      %
      % Usage
      %   fray = ray.focus(location)
      %
      % Parameters
      %   - location (3x1 numeric) -- Location to focus rays towards.

      rdir = cross(ray, location - ray.origin);
      rmat = ott.utils.rotation_matrix(rdir.direction);
      ray = ray.rotate(rmat);
    end
  end

  methods (Hidden)
    function p = getBeamPower(beam)
      % Returns finite power of the ray set
      p = sum(beam.field(:).^2);
    end

    function E = efieldInternal(beam, xyz, varargin)
      % Calculate the E-field
      %
      % Usage
      %   E = beam.efieldInternal(xyz, ...)
      %
      % Optional named arguments
      %   - method (enum) -- Method to use when drawing rays
      %     Can be one of
      %     - inf -- Infinite extend plane waves
      %     - ray_invr -- Power scaled with 1/distance from centre
      %     - invr -- Power scaled with 1/distance from origin
      %     Default: ``'ray_invr'``.

      % Change default parameters
      p = inputParser;
      p.addParameter('method', 'ray_invr');
      p.parse(varargin{:});

      E = efieldInternal@ott.beam.PlaneWave(beam, xyz, ...
          'method', p.results.method);
    end
  end
end

