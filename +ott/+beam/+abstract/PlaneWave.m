classdef PlaneWave < ott.beam.properties.PlaneWave ...
    & ott.beam.abstract.Beam ...
    & ott.beam.utils.InfinitePower
% Abstract representation of a plane wave beam.
% Inherits from :class:`Beam` and :class:`ott.beam.properties.PlaneWave`.
%
% This class creates a single abstract Plane Wave beam.  For more
% efficient arrays of plane wave beams, use :class;`ott.beam.PlaneWave`
% directly.
%
% Supported casts
%   - Beam            -- Default Beam cast, uses PlaneWave
%   - PlaneWave
%   - Ray
%   - vswf.Bsc        -- Default Bsc cast, uses vswf.PlaneWave
%   - vswf.PlaneWave
%   - vswf.PlaneBasis

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = PlaneWave(varargin)
      % Construct a new abstract plane wave representation
      %
      % Usage
      %   beam = PlaneWave(...)
      %
      % Optional named arguments
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

      % Parse parameters
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('direction', []);
      p.addParameter('polarisation', []);
      p.addParameter('origin', []);
      p.addParameter('field', []);
      p.addParameter('vector', []);
      p.addParameter('array_type', 'coherent');
      p.addParameter('like', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Handle defaults from like
      default_direction = [0;0;1];
      default_origin = [0;0;0];
      default_field = 1.0;
      default_polarisation = [1;0;0];
      if ~isempty(p.Results.like)
        if isa(p.Results.like, 'ott.beam.abstract.PlaneWave')
          default_direction = p.Results.like.direction;
          default_origin = p.Results.like.origin;
          default_field = p.Results.like.field;
          default_polarisation = p.Results.like.polarisation;
        end
      end

      % Get which parameters are using defaults (i.e. unset)
      udef_origin = any(strcmpi('origin', p.UsingDefaults));
      udef_direction = any(strcmpi('direction', p.UsingDefaults));
      udef_field = any(strcmpi('field', p.UsingDefaults));
      udef_polarisation = any(strcmpi('polarisation', p.UsingDefaults));
      udef_vector = any(strcmpi('vector', p.UsingDefaults));

      % Get values from origin/direction or vector
      assert(udef_vector || (udef_direction && udef_origin), ...
        'vector parameter incompatible with direction/origin');
      if ~udef_vector
        origin = p.Results.vector.origin;
        direction = p.Results.vector.direction;
      else
        if udef_origin
          origin = default_origin;
        else
          origin = p.Results.origin;
        end
        if udef_direction
          direction = default_direction;
        else
          direction = p.Results.direction;
        end
      end

      % Get Vector to store most
      beam = beam@ott.beam.utils.ArrayType('array_type', p.Results.array_type);
      beam = beam@ott.utils.Vector(origin, direction);
      beam = beam@ott.beam.abstract.Beam(unmatched{:}, ...
          'like', p.Results.like);

      % Store remaining parameters
      if udef_field
        beam.field = default_field;
      else
        beam.field = p.Results.field;
      end
      if udef_polarisation
        beam.polarisation = default_polarisation;
      else
        beam.polarisation = p.Results.polarisation;
      end
    end

    function ray = ott.beam.Ray(plane)
      % Type conversion from plane wave to ray
      ray = ott.beam.Ray('origin', plane.origin, ...
        'polarisation', plane.polarisation, ...
        'field', plane.field, ...
        'direction', plane.direction);
    end

    function bsc = ott.beam.vswf.Bsc(plane, varargin)
      % Convert from a PlaneWave to a Bsc
      %
      % Usage
      %   bsc = ott.beam.vswf.Bsc(planewave, ...)
      %
      % Optional named arguments
      %   - suggestedNmax (numeric) -- Suggested Nmax for the Bsc.
      %     For plane waves, this is the Nmax used for the beam.

      p = inputParser;
      p.addParameter('suggestedNmax', []);
      p.parse(varargin{:});

      % Calculate wave directions
      rtp = ott.utils.xyz2rtp(plane.direction);
      theta = rtp(2, :);
      phi = rtp(3, :);

      % Construct BSC
      bsc = ott.beam.vswf.Plane(theta, phi, 'Nmax', p.Results.suggestedNmax);
    end
  end
end
