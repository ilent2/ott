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
%   - Beam        -- Creates an empty Array
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

    function beam = ott.beam.Beam(beam, varargin)
      beam = castArrayHelper(@ott.beam.Array.like, beam, ...
          'array_type', 'coherent', varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(beam, varargin)
      beam = castArrayHelper(@ott.beam.vswf.Bsc.like, beam, ...
          'array_type', 'coherent', varargin{:});
    end

    function beam = ott.beam.PlaneWave(beam, varargin)
      beam = castArrayHelper(@ott.beam.PlaneWave.like, beam, ...
          'array_type', 'coherent', varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)
      beam = castArrayHelper(@ott.beam.Ray.like, beam, ...
          'array_type', 'coherent', varargin{:});
    end

    function beam = ott.beam.Dipole(beam, varargin)
      beam = castArrayHelper(@ott.beam.Dipole.like, beam, varargin{:});
    end
  end

  methods (Access=protected)
    function beam = castArrayHelper(cast, beam, varargin)
      % Helper for casts

      assert(isa(beam, 'ott.beam.abstract.Empty'), ...
          'First argument must be a abstract.Empty');
      ott.utils.nargoutCheck(beam, nargout);

      beam = cast(beam(1), varargin{:});
    end
  end
end
