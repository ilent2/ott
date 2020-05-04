classdef (Abstract) Abstract < ott.beam.BeamProperties
% Base class for beam descriptions (no fields)
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
    function beam = Abstract(varargin)
      % Construct an abstract beam instance
      %
      % beam = AbstractBeam(...)
      %
      % For optional parameters, see :class:`BeamProperties`.

      beam = beam@ott.beam.BeamProperties(varargin{:});
    end
  end
end

