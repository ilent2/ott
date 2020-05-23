function args = addDefaultParameter(name, val, args)
% Add a default parameter to a cell array of function arguments
%
% Usage
%   args = addDefaultParameter(name, val, varargin)
%
% Parameters
%   - name (char) -- Parameter name
%   - val (any) -- Parameter default value
%   - varargin (cell) -- Cell array of parameters
%
% Returns a cell array of parameters.

  if ~any(strcmpi(name, args))
    args = [args, {name}, {val}];
  end
end
