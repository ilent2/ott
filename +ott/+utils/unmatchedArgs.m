function unmatched = unmatchedArgs(p)
% Gets a cell array of the unmatched arguments in an input parser
%
% Usage
%   unmatched = unmatchedArgs(p)
%   To use the result: ``another_function(unmatched{:});``.

unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];