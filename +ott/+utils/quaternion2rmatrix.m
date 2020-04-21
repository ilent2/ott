function rmatrix = quaternion2rmatrix(quat)
% Convert from quaternion to rotation matrix representation
%
% Usage
%   rmatrix = rmatrix2quaternion(rmatrix)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  assert(size(quat, 1) == 4, 'Quaternion must be 4xN matrix');

  % Assign simple names
  q1 = quat(1, :);
  q2 = quat(2, :);
  q3 = quat(3, :);
  q4 = quat(4, :);

  % Allocate memory
  rmatrix = zeros(3, 3, size(q4, 2));

  rmatrix(1, 1, :) = q1.^2 - q2.^2 - q3.^2 + q4.^2;
  rmatrix(2, 2, :) = -q1.^2 + q2.^2 - q3.^2 + q4.^2;
  rmatrix(3, 3, :) = -q1.^2 - q2.^2 + q3.^2 + q4.^2;

  rmatrix(2, 1, :) = 2.*(q1.*q2 - q3.*q4);
  rmatrix(3, 1, :) = 2.*(q1.*q3 + q2.*q4);
  rmatrix(3, 2, :) = 2.*(q2.*q3 - q1.*q4);
  rmatrix(1, 2, :) = 2.*(q1.*q2 + q3.*q4);
  rmatrix(1, 3, :) = 2.*(q1.*q3 - q2.*q4);
  rmatrix(2, 3, :) = 2.*(q2.*q3 + q1.*q4);

  rmatrix = reshape(rmatrix, 3, []);

end

