function insertGridRow(grid, rowheight, idx, numrows)
% Insert additional grid rows and shuffle everything else down
%
% insertGridRow(grid, rowheight, idx, numrows)
%
%  Inserts numrows new rows with rowheight at index idx into grid.

% Copyright IST Austria 2022
% Written by Isaac Lenton


  % Check value idx
  assert(isnumeric(idx) && idx > 0 && idx <= numel(grid.RowHeight), ...
    'idx must be a valid row index in the grid');
  assert(isnumeric(numrows) && numrows > 0, ...
    'numrows must be positive numeric integer');
  
  % Extend the grid
  grid.RowHeight(:, idx+numrows:end+numrows) = grid.RowHeight(:, idx:end);
  grid.RowHeight(:, idx:idx+numrows-1) = repmat({rowheight}, 1, numrows);
    
  % For each child in the insert row or above, add numrows to the row idx
  for ii = 1:numel(grid.Children)
    if grid.Children(ii).Layout.Row >= idx
      grid.Children(ii).Layout.Row = grid.Children(ii).Layout.Row + numrows;
    end
  end
end
