classdef Parity
% Declares a parity beam property
%
% Properties
%   - parity      -- Parity of beam (either 'even' or 'odd')

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    parity       % Parity of beam (either 'even' or 'odd')
  end

  properties (Hidden, SetAccess=protected)
    parityInternal
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.parity(beam, val)
      assert(sum(strcmpi(val, {'even', 'odd'})) == 1, ...
          'parity must be ''even'' or ''odd''');
      beam.parityInternal = val;
      beam.data = [];
    end
    function val = get.parity(beam)
      val = beam.parityInternal;
    end
  end
end
