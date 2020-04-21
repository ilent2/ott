function quat = rmatrix2quaternion(rmatrix)
% Convert from rotation matrix to quaternion rotation representation
%
% Usage
%   quat = rmatrix2quaternion(rmatrix)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  if size(rmatrix, 2) > 3
    rmatrix = reshape(rmatrix, 3, 3, []);
  end

  q4 = 0.5*sqrt(1 + rmatrix(1, 1, :) + rmatrix(2, 2, :) + rmatrix(3, 3, :));

  % Check valid values
  if any(q4 == 0)
    warning('ott:utils:rmatrix2quaternion:zero_division', ...
        'Could not convert some angles');
  end

  q1 = 1./(4.*q4) .* (rmatrix(2, 3, :) - rmatrix(3, 2, :));
  q2 = 1./(4.*q4) .* (rmatrix(3, 1, :) - rmatrix(1, 3, :));
  q3 = 1./(4.*q4) .* (rmatrix(1, 2, :) - rmatrix(2, 1, :));

  % Package output
  quat = [q1; q2; q3; q4];
  quat = squeeze(quat);

end
