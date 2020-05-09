classdef Empty < ott.beam.abstract.Beam & ott.beam.utils.ArrayType
% An empty beam used for empty/default array elements in beam arrays.
% Inherits from :class:`Beam` and :class:`ott.beam.utils.ArrayType`.
%
% When this array is combined with another array beam it dissolves.
% For example, `[Empty, PlaneWave]` would produce `PlaneWave`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Should this move to ott.beam.utils?
% TODO: Should this inherit from ott.beam.abstract.Beam?  Why?

  methods (Hidden)
    function P = getBeamPower(beam)
      % get method called by dependent property power
      P = 0.0;
    end

    function sz = size(beam, dim)
      % Get array size (empty)
      if nargin == 2
        sz = 0;
      else
        sz = [0, 0];
      end
    end

    function varargout = catInternal(beam, varargin)
      error('Empty beams cant be concatenated');
    end

    function varargout = plusInternal(beam, varargin)
      error('Empty beams cant be concatenated');
    end

    function varargout = subsrefInternal(beam, varargin)
      error('Empty beams have no content');
    end

    function varargout = subsasgnInternal(beam, varargin)
      error('Empty beams have no content');
    end
  end
end
