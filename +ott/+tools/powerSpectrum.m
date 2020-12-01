function h = powerSpectrum(varargin)
% Calculate a power spectrum for the given time series data
%
% Currently requires number of points to be even evenly spaced.
% This may change in future.
%
% Usage
%   h = powerSpectrum(t, x, ...)
%
%   h = powerSpectrum(axis, ...)
%   As above, but also specify the plot axis to use.
%
% Parameters
%   - t (N numeric) -- Time for each point.  Must be evenly spaced for now.
%   - x (MxN numeric) -- Position for each point.  Generates a power
%     spectrum for each row in x.
%
% Returns
%   - h (optional) -- Plot handles for generated lines.
%
% Optional named arguments
%   - frequency_units (char) -- Units to use for x label.
%     Default: ``'Hz'``.
%
%   - power_units (char) -- Units to use for y label.
%     Default: ``'a.u.'`` (may change in future).

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  % Get axis handle
  if isa(varargin{1}, 'matlab.graphics.axis.Axes')
    oax = varargin{1};
    varargin = varargin{2:end};
  else
    oax = gca();
  end

  % Parse remaining inputs
  p = inputParser;
  p.addRequired('t');
  p.addRequired('x');
  p.addParameter('frequency_units', 'Hz');
  p.addParameter('power_units', 'a.u.');
  p.parse(varargin{:});
  
  x = p.Results.x;
  t = p.Results.t;
  
  % Ensure even
  if mod(numel(t), 2) ~= 0
    warning('Reducing number of points to nearest even number');
    x = x(:, 1:end-1);
    t = t(1:end-1);
  end

  xw = fft(x, [], 2);
  xw = xw(:, 1:end/2+1);

  Nelm = numel(t);
  Fs = 1./diff(t(1:2));
  f = 0:Fs/Nelm:Fs/2;

  h = loglog(oax, f, abs(xw.^2));
  xlabel(oax, ['Frequency [', p.Results.frequency_units ']']);
  ylabel(oax, ['Power [', p.Results.power_units ']']);

end
