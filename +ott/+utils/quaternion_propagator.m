function qprop = quaternion_propagator(rotation_vector)
% Calculates the quaternion propagation matrix from a rotation vector
%
% Calculates :math:`\frac{dq}{dt}` which can be used to find::
%
%   q^{t+1} = q^{t} + \frac{dq}{dt} \Delta_t
%
% Usage
%   qprop = quaternion_propagator(rotation_vector)
%   Calculates the 4x4 propagation matrix from the 3 element rotation vector.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  rt = rotation_vector(3);
  qt = rotation_vector(2);
  pt = rotation_vector(1);

  qprop = [0, rt, -qt, pt;
           -rt, 0, pt, qt;
           qt, -pt, 0, rt;
           -pt, -qt, -rt, 0];

end

