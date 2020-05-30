classdef Incoherent < ott.beam.properties.AbstractArray ...
    & ott.beam.abstract.Beam
% Incoherent abstract beam arrays.
% Inherits from :class:`ott.beam.properties.AbstractArray` and :class:`Beam`.
%
% Methods
%   - contains          -- Query if array contains specific array_type
%
% Supported casts
%   - Beam              -- Default Beam cast, uses Incoherent
%   - Incoherent
%   - Coherent          -- Not recommended, raises a warning
%   - Array             -- Not recommended, raises a warning
%   - abstract.Coherent -- Convert to incoherent array
%   - abstract.Array    -- Convert to generic array
%   - Bsc
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

      if strcmpi(array_type, 'incoherent')
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
      beam.array_type = 'incoherent';
    end

    function beam = ott.beam.PlaneWave(beam, varargin)
      % Cast to PlaneWave
      beam = castHelper(@ott.beam.PlaneWave.like, ...
          beam, varargin{:});
      beam.array_type = 'incoherent';
    end

    function beam = ott.beam.Ray(beam, varargin)
      % Cast to Ray
      beam = castHelper(@ott.beam.Ray.like, ...
          beam, varargin{:});
      beam.array_type = 'incoherent';
    end
  end
end
