function [velocity, xedges, yedges, counts] = velocityField(t, xy, varargin)
% Calculate 2-D velocity field from position and time data
%
% Usage
%   [velocity, xedges, yedges, counts] = velocityField(t, x, ...)
%
% Returns
%   - velocity (AxBx2 numeric) -- Average velocity at each histogram location.
%   - xedges (A numeric) -- Histogram bin edges.
%   - yedges (B numeric) -- Histogram bin edges.
%   - counts (AxB numeric) -- Histogram counts.
%
% Parameters
%   - t (N numeric) -- Time data.
%
%   - x (2xN numeric) -- Cartesian positions ``[x; y]``.
%
% Unmatched parameters are passed to ``histcounts2``.  For example, to
% specify bin edges, use::
%
%   [velocity, ...] = velocityField(t, x, XEDGES, YEDGES)

  % Calculate positions and velocities with central differences
  avgPos = movmean(xy, 2, 2, 'Endpoints', 'discard');
  vels = diff(xy, 1, 2)./diff(t, 1, 2);

  % Generate histogram data
  [counts, xedges, yedges, xbins, ybins] = ...
      histcounts2(avgPos(1, :), avgPos(2, :), varargin{:});


  % Calculate mean velocities
  velocity = zeros([size(counts), 2]);
  for ii = 1:numel(yedges)-1
    for jj = 1:numel(xedges)-1
      velocity(jj, ii, :) = mean(vels(:, ybins == ii & xbins == jj), 2);
    end
  end

end

