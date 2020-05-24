classdef Empty < ott.beam.properties.Empty ...
    & ott.beam.abstract.Beam
% Empty abstract beam.
% Inherits from :class:`ott.beam.properties.Empty` and :class:`Beam`.
%
% An empty beam describes a optical medium, position and rotation
% but has no power and no fields.  Empty beams are the default scalar
% element for abstract beam hetrogeneous arrays.
%
% Supported casts
%   - Beam        -- Default beam cast, uses Empty
%   - Empty
%   - Array       -- (Inherited from abstract.Beam)
%   - Coherent    -- (Inherited from abstract.Beam)
%   - Incoherent  -- (Inherited from abstract.Beam)
%   - PlaneWave   -- Plane-wave array with no beams
%   - Bsc         -- Bsc array with no beams
%   - Ray         -- Empty array with no beams
%   - Dipole      -- Empty array with no beams

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Empty(varargin)
      % Construct an empty beam.
      %
      % Empty beams still have beam properties but have zero power and
      % no fields.  For properties, see :class:`ott.beam.properties.Empty`.

      beam = beam@ott.beam.properties.Empty(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      beam = ott.beam.Empty(varargin{:});
    end

    function beam = ott.beam.Empty(beam, varargin)
      % Cast to a Empty or coherent array of empties

      assert(isa(beam, 'ott.beam.abstract.Empty'), ...
          'First argument must be a abstract.Empty');

      if numel(beam) > 1
        beam = ott.beam.Coherent(beam, varargin{:})
      else
        beam = ott.beam.Empty.like(beam, varargin{:});
      end
    end

    function beam = ott.beam.vswf.Bsc(beam, varargin)
      % Cast to an empty Bsc, discard meta for arrays of beams

      assert(isa(beam, 'ott.beam.abstract.Empty'), ...
          'First argument must be a abstract.Empty');

      beam = ott.beam.vswf.Bsc.like(beam(1), varargin{:});
    end

    function beam = ott.beam.PlaneWave(beam, varargin)
      % Cast to an empty PlaneWave, discard meta for arrays of beams

      assert(isa(beam, 'ott.beam.abstract.Empty'), ...
          'First argument must be a abstract.Empty');

      beam = ott.beam.PlaneWave.like(beam(1), varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Cast to an empty PlaneWave, discard meta for arrays of beams

      assert(isa(beam, 'ott.beam.abstract.Empty'), ...
          'First argument must be a abstract.Empty');

      beam = ott.beam.Ray.like(beam(1), varargin{:});
    end

    function beam = ott.beam.Dipole(beam, varargin)
      % Cast to an empty Dipole, discard meta for arrays of beams

      assert(isa(beam, 'ott.beam.abstract.Empty'), ...
          'First argument must be a abstract.Empty');

      beam = ott.beam.Dipole.like(beam(1), varargin{:});
    end
  end
end
