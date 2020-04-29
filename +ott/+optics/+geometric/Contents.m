% ott.optics.geometric Classes for geometric optics calculations
%
% TODO

% Ray -- Object describing geometric rays
%   Supports calculating scattering by a plane
%   Supports calculating momentum?
%   Provides a scatter method (taking normals as inputs)
%
% Beam?
%   In OTGO these described the paraxial beams
%   We just need a way to construct Ray objects for different beams
%   Maybe we should provide a .simple constructor to Ray and
%   wrap a beam object described elsewhere?
%
% Particle -- Object describing a geometric optics particle
%   Wraps a shape object
%   Provides a force and torque method and a scatter method
%   To be consisten with bsc, should we only have scatter as part of Ray?
%
% System -- For use with dynamics systems
%   Wraps a Ray and Particle object
%   Provides the force(position, rotation, ...) method needed for dynamics
%   Need to review how other optics method define/will define a system
%
% We really need to review other optics packages in order to provide
% a consistent feel to the toolbox.
% Where are we going to put Davis beams and paraxial beam representations?

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

