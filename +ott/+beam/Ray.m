classdef Ray < ott.beam.properties.PlaneWaveArray ...
    & ott.beam.properties.AnyArrayType ...
    & ott.beam.Beam
% Specialisation of a PlaneWave for Ray-optics beams
%
% The main differences between :class:`PlaneWave` and :class:`Ray` are
% finite power and more ray-orientated defaults for methods.
%
% Static methods
%   - like            -- Create a beam like another
%   - FromDirection   -- Construct from direction/polarisation
%   - DirectionSet    -- Construct a directionSet
%   - empty           -- Construct an empty array
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

  properties (Dependent)
    power       % Power calculated from the rays
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.properties.AnyArrayType.likeProperties(other, args);
      args = ott.beam.properties.PlaneWaveArray.likeProperties(other, args);
      args = ott.beam.Beam.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Gaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.Ray.likeProperties(other, varargin);
      beam = ott.beam.Ray(args{:});
    end

    function beam = FromDirection(varargin)
      % Construct beam from direction/polarisation vectors.
      %
      % Usage
      %   beam = FromDirection(origin, direction, polarisation, field)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - origin (3xN numeric) -- Origin (for phase offset) of wave.
      %   - direction (3xN numeric) -- Propagation direction of wave.
      %   - polarisation (3xN numeric) -- Primary polarisation direction.
      %   - field (2xN numeric) -- Field in two polarisation directions.
      %
      % Additional parameters passed to base.

      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('direction', [], @isnumeric);
      p.addOptional('polarisation1', [], @isnumeric);
      p.addOptional('field', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      % Construct direction set
      directionSet = ott.beam.properties.PlaneWave.DirectionSet(...
          p.Reults.direction, p.Results.polarisation1);

      % Construct beam
      beam = ott.beam.Ray(...
          'origin', p.Results.origin, ...
          'directionSet', directionSet, ...
          'field', p.Results.field, ...
          unmatched{:});
    end

    function beam = empty(varargin)
      % Construct an emtpy beam array
      %
      % Usage
      %   beam = ott.beam.Ray.empty(...)
      %
      % Additional parameters are passed to the constructor.

      beam = ott.beam.Ray('directionSet', zeros(3, 0), ...
        'field', zeros(2, 0), 'origin', zeros(3, 0), varargin{:});
    end
  end

  methods
    function bm = Ray(varargin)
      % Construct a collection of Rays
      %
      % Usage
      %   beam = Ray(origin, directionSet, ...)
      %
      % See also :meth:`FromDirection` for construction from
      % direction/polarisation vectors.
      %
      % Parameters
      %   - origin (3xN numeric) -- Plane wave origins.
      %
      %   - directionSet (3x3N numeric) -- Array formed by combining
      %     direction/polarisation vectors into rotation matrices.  The
      %     direction vector should be the last column of the matrix.
      %     Default: ``eye(3)``.
      %
      % Optional named arguments
      %   - field (2xN numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %     Default: ``[1; 1i]``.
      %
      %   - array_type (enum) -- Beam array type.  Can be
      %     'coherent', 'incoherent' or 'array'.  Default: ``'coherent'``.

      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('directionSet', [], @isnumeric);
      p.addParameter('field', [1;1i]);
      p.addParameter('array_type', 'coherent');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bm = bm@ott.beam.properties.AnyArrayType(p.Results.array_type);
      bm = bm@ott.beam.properties.PlaneWaveArray(...
          'origin', p.Results.origin, ...
          'directionSet', p.Results.directionSet, ...
          'field', p.Results.field, ...
          unmatched{:});
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

    % TODO: Visualise far-field methods with points

    function varargout = visualise(vec, varargin)
      % Plots the vector set in 3-D.
      %
      % Uses the quiver function to generate a visualisation of the
      % vector set.
      %
      % Usage
      %   h = beam.visualise(...)
      %
      % Optional named arguments
      %   - Scale (numeric) -- rescales the coordinates and components
      %     of the vector before plotting.  Can either be a scalar
      %     or vector ``[S1, S2]`` specifying separate scaling for the
      %     coordinates and components.  Default: ``[1, 1]``.
      %
      %   - ray_lengths (1xN numeric) -- Array of ray lengths for
      %     plotting finite length rays.
      %
      % Any unmatched named arguments are applied to the plot handles
      % returned by the quiver function calls.

      % Parse inputs
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('Scale', [1, 1]);
      p.addParameter('ray_lengths', []);
      p.addParameter('show_polarisation', true);
      p.parse(varargin{:});

      S1 = p.Results.Scale(1);
      S2 = p.Results.Scale(2);
      S2pol = p.Results.Scale(2);

      if ~isempty(p.Results.ray_lengths)
        S2 = p.Results.ray_lengths;
        assert(size(S2, 1) == 1, 'ray_lengths must be 1xN');
      end

      isholdon = ishold();

      % Generate plot of directions
      h = quiver3(S1*vec.origin(1, :), S1*vec.origin(2, :), ...
          S1*vec.origin(3, :), S2.*vec.direction(1, :), ...
          S2.*vec.direction(2, :), S2.*vec.direction(3, :), 0);

      if p.Results.show_polarisation

        if ~isholdon
          hold('on');
        end

        % Generate plot of polarisations
        h(2) = quiver3(S1*vec.origin(1, :), S1*vec.origin(2, :), ...
            S1*vec.origin(3, :), S2pol.*vec.polarisation1(1, :), ...
            S2pol.*vec.polarisation1(2, :), S2pol.*vec.polarisation1(3, :), 0);

        if ~isholdon
          hold('off');
        end

      end

      % Apply unmatched arguments to plot handle
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];
      if ~isempty(unmatched)
        set(h, unmatched{:});
      end

      % Assign outputs
      if nargout > 0
        varargout{1} = h;
      end
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, xyz, varargin)
      % Cast to PlaneWave for visualisation

      % Change default parameters
      p = inputParser;
      p.addParameter('method', 'ray_invr');
      p.parse(varargin{:});

      beam = ott.beam.PlaneWave(beam, 'position', [0;0;0]);
      E = beam.efield(xyz, 'method', p.Results.method);
    end

    function H = hfieldInternal(beam, xyz, varargin)
      % Cast to PlaneWave for visualisation

      % Change default parameters
      p = inputParser;
      p.addParameter('method', 'ray_invr');
      p.parse(varargin{:});

      beam = ott.beam.PlaneWave(beam, 'position', [0;0;0]);
      H = beam.hfield(xyz, 'method', p.Results.method);
    end

    function E = efarfieldInternal(beam, rtp, varargin)
      % Cast to PlaneWave for visualisation

      beam = ott.beam.PlaneWave(beam, 'position', [0;0;0]);
      E = beam.efarfield(rtp, 'method', p.Results.method);
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Cast to PlaneWave for visualisation

      beam = ott.beam.PlaneWave(beam, 'position', [0;0;0]);
      H = beam.hfarfield(rtp, 'method', p.Results.method);
    end
  end

  methods % Getters/setters
    function p = get.power(beam)
      % Returns finite power of the ray set
      p = beam.intensity;
    end
  end
end

