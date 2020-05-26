classdef Coherent < ott.beam.properties.CoherentArrayType
% Pseudonym class for coherent arrays.
% Inherits from :class:`properties.CoherentArrayType`.
%
% This class acts to identify coherent arrays and enable casts directly
% to coherent arrays.  To construct a new coherent array, use the Array
% class instead.
%
% Usage
%
%   coherent_array = ott.beam.Coherent(other_array)
%   Casts the array type to a coherent array.  See cast method in
%   `other_array` for details about implementation/usage.
%
%   b = isa(array_type, 'ott.beam.Coherent')
%   Returns true when the 'array_type' property is 'coherent'.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Access=private)
    function beam = Coherent(varargin)
      % Private constructor, this type isn't intended to be used directly.
      % Use the Array class instead.
    end
  end
end
