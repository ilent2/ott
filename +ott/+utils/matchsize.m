function varargout = matchsize(varargin)
% Checks that all vector inputs have the same number of rows.
%
% Usage
%   [A,B,...] = matchsize(A,B,...) checks inputs have same number of rows,
%   and expands single-row inputs by repetition to match the input row number.
%
% Parameters
%   - A,B,... (numeric)  -- Numeric arrays whose number of rows
%     are to be matched.
%
% Example
%   The following example shows has two inputs, a scalar and a row vector.
%   The scalar is duplicated to match the length of the row vector::
%
%     A = 5;
%     B = [1; 2; 3];
%
%     [A,B] = matchsize(A, B);
%     disp(A)  % -> [5; 5; 5]
%     disp(B)  % -> [1; 2; 3]

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
