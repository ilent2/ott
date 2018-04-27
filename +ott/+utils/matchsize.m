function varargout = matchsize(varargin)
% MATCHSIZE checks that all vector inputs have the same number of rows
%
% [A,B,...] = MATCHSIZE(A,B,...) checks inputs have same number of rows,
% and expands single-row inputs by repetition to match the input row number.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Loop over inputs to determine maximum size
nmax = 0;
for ii = 1:length(varargin)
  nmax = max(nmax, size(varargin{ii}, 1));
end

% Loop over inputs to check size/generate output
for ii = 1:length(varargin)
  nrows = size(varargin{ii}, 1);
  if nrows == 1
    varargout{ii} = repmat(varargin{ii}, nmax, 1);
  elseif nrows ~= nmax
    error('ott:matchsize:size_mismatch', ...
        'Number of rows in inputs must be one or equal.');
  else
    varargout{ii} = varargin{ii};
  end
end
