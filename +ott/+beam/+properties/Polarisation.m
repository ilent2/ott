classdef Polarisation
% Declares a polarisation basis property.
%
% Properties
%   - polbasis    -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield    -- (2 numeric) Field in the theta/phi or x/y directions

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    polbasis    % (enum) Polarisation basis ('polar' or 'cartesian')
    polfield    % (2 numeric) Field in the theta/phi or x/y directions
  end

  properties (Hidden, SetAccess=protected)
    polbasisInternal
    polfieldInternal
  end

  properties (Abstract)
    data
  end

  methods % Getters/setters
    function beam = set.polbasis(beam, val)
      assert(sum(strcmpi(val, {'polar', 'cartesian'})) == 1, ...
          'polbasis must be ''polar'' or ''cartesian''');
      beam.polbasisInternal = val;
      beam.data = [];
    end
    function val = get.polbasis(beam)
      val = beam.polbasisInternal;
    end

    function beam = set.polfield(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'polfield must be 2 element numeric vector');
      beam.polfieldInternal = [val(1), val(2)];
      beam.data = [];
    end
    function val = get.polfield(beam)
      val = beam.polfieldInternal;
    end
  end
end
