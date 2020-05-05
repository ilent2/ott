classdef (Abstract) Beam < ott.beam.Properties
% Base class for beam descriptions (no fields).
% Inherits from :class:`ott.beam.Properties`.
%
% Classes that inherit from this class include a description of the
% beam but do not necessarily include methods to calculate the fields.
%
% Abstract methods
%   - getBeamPower      -- get method called by dependent property power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Beam(varargin)
      % Construct an abstract beam instance
      %
      % beam = Beam(...)
      %
      % For optional parameters, see :class:`Properties`.

      beam = beam@ott.beam.Properties(varargin{:});
    end
  end
end

