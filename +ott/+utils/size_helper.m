function varargout = size_helper(sz, dim)
% Helper for implementing matlab size-like functionality
%
% Usage
%   varargout = size_helper(size(quantity))

if nargin == 2
  varargout{1} = sz(dim);
else
  if nargout > 1
    if nargout > numel(sz)
      sz = [sz, ones(1, nargout - numel(sz))];
    elseif nargout < numel(sz)
      sz = [sz(1:nargout-1), prod(sz(nargout:end))];
    end

    [varargout{:}] = deal(sz);
  else
    varargout{1} = sz;
  end
end
