function onecol = col3to1(threecol)
% Reshape a 3 column matrix into a 1 column vector
%
% onecol = col3to1(threecol) reshape the threecol to a onecol:
%   threecol = [x1, y1, z1; x2, y2, z2; ... ; xN, yN, zN]
%   onecol = [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN]

onecol = reshape(threecol.', [numel(threecol), 1]);

