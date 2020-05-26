classdef Coherent < ott.beam.properties.AbstractArray ...
    & ott.beam.abstract.Beam
% Coherent abstract beam arrays.
% Inherits from :class:`ott.beam.properties.AbstractArray` and :class:`Beam`.
%
% Coherent abstract beam arrays can be created with the array syntax,
% for example ``[beam1, beam2, beam3]``.  Therefore, this class is
% somewhat redundant for creating coherent arrays, it is included for
% consistently and creating arrays of coherent sub-arrays.
%
% Methods
%   - contains        -- Query if array contains specific array_type
%
% Supported casts
%   - Beam                  -- Default Beam cast, uses Coherent
%   - Coherent
%   - Incoherent            -- Not recommended, raises a warning
%   - Array                 -- Not recommended, raises a warning
%   - abstract.Incoherent   -- Convert to incoherent array
%   - abstract.Array        -- Convert to generic array
%   - Bsc
%   - Dipole
%   - PlaneWave
%   - Ray

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function b = contains(beam, array_type)
      % Query if array contains specific array_type
      %
      % Usage
      %   b = beam.contains(array_type)
      %
      % Parameters
      %   - array_type (enum) -- An array type, must be one of
      %     'array', 'coherent' or 'incoherent'.

      if strcmpi(array_type, 'coherent')
        b = true;
      elseif numel(beam) > 1
        b = contains@ott.beam.abstract.Beam(beam, array_type);
      else
        b = false;
        for ii = 1:numel(beam.beams)
          b = b | beam.beams(ii).contains(array_type);
        end
      end
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to Incoherent
      beam = ott.beam.Incoherent(varargin{:});
    end

    function beam = ott.beam.Incoherent(beam, varargin)
      % Cast to Incoherent
      beam = castHelper(@ott.beam.Incoherent.like, beam, varargin{:});
    end

    function beam = ott.beam.Coherent(beam, varargin)
      % Cast to Incoherent
      beam = castHelper(@ott.beam.Coherent.like, beam, varargin{:});
      warning('ott:beam:abstract:Incoherent:different_type', ...
          'Cast to different array type, perhaps you meant Incoherent');
    end

    function beam = ott.beam.Array(beam, varargin)
      % Cast to Incoherent
      beam = castHelper(@ott.beam.Array.like, beam, varargin{:});
      warning('ott:beam:abstract:Incoherent:different_type', ...
          'Cast to different array type, perhaps you meant Incoherent');
    end

    function beam = ott.beam.abstract.Coherent(beam, varargin)
      % Cast to abstract.Coherent
      beam = castHelper(@ott.beam.abstract.Coherent.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.abstract.Array(beam, varargin)
      % Cast to abstract.Array
      beam = castHelper(@ott.beam.abstract.Array.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(beam, varargin)
      % Cast to vswf.Bsc
      beam = castHelper(@ott.beam.vswf.Bsc.like, ...
          beam, varargin{:});
      beam.array_type = 'coherent';
    end

    function beam = ott.beam.Dipole(beam, varargin)
      % Cast to Dipole
      beam = castHelper(@ott.beam.Dipole.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.PlaneWave(beam, varargin)
      % Cast to PlaneWave
      beam = castHelper(@ott.beam.PlaneWave.like, ...
          beam, varargin{:});
      beam.array_type = 'coherent';
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Cast to Ray
      beam = castHelper(@ott.beam.Ray.like, ...
          beam, varargin{:});
      beam.array_type = 'coherent';
    end
  end

  methods (Access=private)
    function beam = castHelper(cast, beam, varargin)
      % Helper for casts

      assert(isa(beam, 'ott.beam.abstract.Coherent'), ...
          'First argument should be a abstract.Coherent');

      if numel(beam) > 1
        oldbeam = beam;
        beam = ott.beam.Array('coherent', size(beam));
        for ii = 1:numel(oldbeam)
          beam(ii) = cast(beam, varargin{:});
        end
      else
        beam = cast(beam, varargin{:});
      end
    end
  end
end
