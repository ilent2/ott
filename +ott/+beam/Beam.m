classdef (Abstract) Beam < ott.beam.utils.BeamInterface
% Abstract base class for optical tweezers toolbox beam implementations.
% Inherits from :class:`ott.beam.utils.BeamInterface`.
%
% This class acts to identify beam implementations.  Beams can be
% created using subclasses of this object or by casting abstract beams
% to the beam object (see :class:`ott.beam.abstract.Beam`).
%
% For a list of function/properties see :class:`ott.beam.utils.BeamInterface`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Beam(varargin)
      % Construct a new beam object
      %
      % Usage
      %   beam = beam@ott.beam.Beam(...)

      beam = beam@ott.beam.utils.BeamInterface(varargin{:});
    end
  end
end
