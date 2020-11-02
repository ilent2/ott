classdef Profile
% Declares a profile property for Annular and Pm beams
%
% Properties
%   - profile       -- Beam profile (function_handle | beam | 'uniform')

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    profile     % Beam profile (function_handle | beam | 'uniform')
  end

  properties (Hidden, SetAccess=protected)
    profileInternal
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.profile(beam, val)
      assert((ischar(val) && strcmpi(val, 'uniform')) ...
          || isa(val, 'ott.beam.Beam') || isa(val, 'ott.bsc.Bsc') ...
          || isa(val, 'function_handle'), ...
          'profile must be ''uniform'', function_handle or beam');
      beam.profileInternal = val;
      beam.data = [];
    end
    function val = get.profile(beam)
      val = beam.profileInternal;
    end
  end
end

