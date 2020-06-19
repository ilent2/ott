classdef StokesData < ott.drag.Stokes
% Stokes drag tensor definition with explicit forward and inverse data
%
% Properties
%   - inverse       -- (6x6 numeric) inverse drag tensor
%   - forward       -- (6x6 numeric) forward drag tensor
%
% Additional methods/properties inherited from :class:`Stokes`.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    inverseInternal
    forwardInternal
  end

  methods
    function drag = StokesData(varargin)
      % Construct new Stokes drag tensor instance with explicit data
      %
      % Usage
      %   drag = StokesDrag(forward, inverse, ...)
      %
      % Parameters
      %   - forward (6x6 numeric) -- Forward drag tensor
      %   - inverse (6x6 numeric) -- Inverse drag tensor
      %
      % If either of `forward` or `inverse` are omitted, they are calculated
      % from the corresponding tensor.
      %
      % Parameters can also be passed as named arguments.
      % Unmatched parameters are passed to :class:`Stokes`.

      p = inputParser;
      p.addOptional('forward', []);
      p.addOptional('inverse', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      drag = drag@ott.drag.Stokes(unmatched{:});

      if isempty(p.Results.forward) && isempty(p.Results.inverse)
        drag.forward = eye(6);
        drag.inverse = eye(6);
      elseif isempty(p.Results.forward)
        drag.inverse = p.Results.inverse;
        drag.forward = inv(drag.inverse);
      elseif isempty(p.Results.inverse)
        drag.forward = p.Results.forward;
        drag.inverse = inv(drag.forward);
      else
        drag.forward = p.Results.forward;
        drag.inverse = p.Results.inverse;
      end
    end
  end

  methods % Getters/setters
    function drag = set.forwardInternal(drag, val)
      assert(ismatrix(val) && all(size(val) == [6, 6]) && isnumeric(val), ...
        'forward must be 6x6 numeric matrix');
      drag.forwardInternal = val;
    end

    function drag = set.inverseInternal(drag, val)
      assert(ismatrix(val) && all(size(val) == [6, 6]) && isnumeric(val), ...
        'inverse must be 6x6 numeric matrix');
      drag.inverseInternal = val;
    end
  end
end

