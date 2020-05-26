classdef Incoherent < ott.beam.properties.IncoherentArrayType
% Pseudonym class for coherent arrays.
% Inherits from :class:`properties.IncoherentArrayType`.
%
% This class acts to identify incoherent arrays and enable casts directly
% to incoherent arrays.  To construct a new incoherent array, use the Array
% class instead.
%
% Usage
%
%   incoherent_array = ott.beam.Incoherent(other_array)
%   Casts the array type to a incoherent array.  See cast method in
%   `other_array` for details about implementation/usage.
%
%   b = isa(array_type, 'ott.beam.Incoherent')
%   Returns true when the 'array_type' property is 'incoherent'.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Access=private)
    function beam = Incoherent(varargin)
      % Private constructor, this type isn't intended to be used directly.
      % Use the Array class instead.
    end
  end
end
