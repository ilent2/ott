function threecol = col1to3(onecol)
% Reshape a 1 column matrix into a 3 column matrix
%
% threecol = col1to3(onecol) reshape the onecol to a threecol:
%   threecol = [x1, y1, z1; x2, y2, z2; ... ; xN, yN, zN]
%   onecol = [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN]

threecol = reshape(onecol, [3, numel(onecol)/3]).';

