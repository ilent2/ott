% ott.dynamics Methods for simulating particle dynamics
%
% Base classes
%   Method          -- Base class for optical tweezers dynamics methods
%   Milstein        -- Base class for Milstein dynamics step size methods
%
% Dynamics methods
%   MatlabOde       -- Solve ODE using one of Matlab's ODE solvers
%   Euler           -- Simple implementation of Euler's method
%   MaxStep         -- Estimates maximum step based on diffusion and force
%   TwoStepMilstein -- Estimates error using two half steps
%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

