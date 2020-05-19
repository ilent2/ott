classdef (Abstract) Beam < ott.beam.Properties
% Base class for beam descriptions (no fields).
% Inherits from :class:`ott.beam.Properties`.
%
% Classes that inherit from this class include a description of the
% beam but do not necessarily include methods to calculate the fields.
%
% Methods
%   - cat               -- Concatenate two beams (creates a Array)
%
% Abstract methods
%   - getBeamPower      -- get method called by dependent property power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Perhaps we need ott.beam.properties and ott.beam.abstract,
%   Abstract beams have properties but are abstract.  Full beams
%   are not abstract but have properties.
%
%   ott.beam.Beam did not (until recently) inherit from
%   ott.beam.abstract.Beam, this meant we could check
%   isa(beam, 'abstract') to determine if the beam had Beam-like properties.
%
%   Abstract beams could probably have their own hetrogeneous
%   array type.  Properties would be shared by Beams and abstract.
%   Abstract is convertible to beams.

  methods
    function beam = Beam(varargin)
      % Construct an abstract beam instance
      %
      % beam = Beam(...)
      %
      % For optional parameters, see :class:`Properties`.

      beam = beam@ott.beam.Properties(varargin{:});
    end

    function beam = horzcat(varargin)
      % Concatenate beam objects.
      %
      % Usage
      %   beam = [beam1, beam2, ...]
      %   Defers to cat(2, ...).

      beam = cat(2, varargin{:});
    end

    function beam = vertcat(varargin)
      % Concatenate beam objects.
      %
      % Usage
      %   beam = [beam1; beam2; ...]
      %   Defers to cat(1, ...).

      beam = cat(1, varargin{:});
    end

    function beam = cat(dim, varargin)
      % Concatenate two beam objects.
      %
      % By default, this creates a new beam Array.
      % If you class has a different behaviour, either overload this
      % function or use the :class:`ott.beam.utils.ArrayType` class.
      %
      % ott.beam.abstract.Empty beams are removed from the array.
      %
      % Usage
      %   beam = cat(dim, beam1, beam2, beam3, ...)

      sz = ones(1, dim+1);
      sz(dim) = numel(varargin);
      beam = catNewArray('array', sz, varargin{:});
    end
  end

  methods (Hidden)
    function beam = catNewArray(array_type, sz, varargin)
      % Construct a new Array for this beam type.
      %
      % This is the default method, which constructs an abstract
      % array of beams.  ott.beam.Beam overloads this method.
      % If your beam type has a custom array type, overload this method.

      beam = ott.beam.abstract.Array(array_type, sz, varargin{:});
    end
  end
end

