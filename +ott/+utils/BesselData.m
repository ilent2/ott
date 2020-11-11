classdef BesselData
% Stores pre-computed Bessel data for repeated field calculations
%
% Properties
%   - un       -- Unique list of n
%   - ur       -- Unique list of r
%   - hn       -- Value
%   - dhn      -- Derivative
%   - func     -- Function use for Bessel calculation
%
% Methods
%   - BesselData    -- Constructor
%   - evaluate      -- Evaluate Bessel function or use stored data
%   - updateN       -- Update stored n values
%   - updateR       -- Update stored r values

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    un        % Unique list of n
    ur        % Unique list of r
    hn        % Value
    dhn       % Derivative
    func      % Function use for Bessel calculation
  end

  methods
    function data = BesselData(func)
      % Construct a new BesselData instance
      %
      % Usage
      %   data = BesselData(func)
      %
      % Parameters
      %   - func (function handle) -- The function to use for calculations.

      data.func = func;
      data.un = [];
      data.ur = [];
      data.hn = [];
      data.dhn = [];
    end

    function [hn, dhn, data] = evaluate(data, n, kr)
      % Evaluate the Bessel function or use stored data
      %
      % Usage
      %   [hn, dhn, data] = data.evaluate(n, kr)
      %
      % Parameters
      %   - n (N numeric) -- mode indices.
      %   - kr (M numeric) -- positions multiplied by wavenumber.
      %
      % Returns
      %   - hn, dhn (NxM numeric) -- Results from evaluating func.
      
      assert(isvector(kr) && isnumeric(kr), ...
          'kr must be numeric vector');
      assert(isvector(n) && isnumeric(n), ...
          'n must be numeric vector');

      % Update data set
      data = data.updateN(n);
      data = data.updateR(kr);

      % Get indices for elements
      [~, dn] = ismember(n, data.un);
      [~, dr] = ismember(kr, data.ur);
      
      keep = dr ~= 0;
      
      % Allocate output
      hn = zeros(numel(n), numel(kr));
      dhn = hn;
      
      % Place nans for invalid values
      hn(:, ~keep) = nan;
      dhn(:, ~keep) = nan;

      % Retrieve stored data
      hn(:, keep) = data.hn(dn, dr(keep));
      dhn(:, keep) = data.dhn(dn, dr(keep));

    end

    function data = updateN(data, n)
      % Update stored N values
      %
      % Usage
      %   data = data.updateN(n)

      ott.utils.nargoutCheck(data, nargout);

      % Find values to update
      n = unique(n);
      [~, dn] = ismember(n, data.un);
      new_n = n(dn == 0);

      % Update values
      N = size(data.hn, 1);
      data.un = [data.un; new_n(:)];
      data.hn = [data.hn; zeros(length(new_n), size(data.hn, 2))];
      data.dhn = [data.dhn; zeros(length(new_n), size(data.hn, 2))];
      for ii = 1:length(new_n)
        [data.hn(N+ii, :), data.dhn(N+ii, :)] = data.func(new_n(ii), data.ur);
      end
    end

    function data = updateR(data, kr)
      % Update stored r values
      %
      % Usage
      %   data = data.updateR(kr)

      ott.utils.nargoutCheck(data, nargout);

      % Find values to update
      kr(isnan(kr)) = [];
      kr = unique(kr);
      [~, dkr] = ismember(kr, data.ur);
      new_kr = kr(dkr == 0);

      % Update values
      N = size(data.hn, 2);
      idx = N + (1:numel(new_kr));
      data.ur = [data.ur; new_kr(:)];
      data.hn = [data.hn, zeros(size(data.hn, 1), numel(new_kr))];
      data.dhn = [data.dhn, zeros(size(data.hn, 1), numel(new_kr))];
      for ii = 1:numel(data.un)
        [data.hn(ii, idx), data.dhn(ii, idx)] = data.func(data.un(ii), new_kr);
      end
    end
  end
end
