% ott.dynamics Methods for simulating particle dynamics
%
% Base classes
%   Method          -- Base class for dynamics methods
%   Milstein        -- Base class for Milstein dynamic step size methods
%   DynamicsSystem  -- Base class for describing dynamics systems
%
% Methods
%   MatlabOde       -- Solve ODE using one of Matlab's ODE solvers
%   Euler           -- Simple implementation of Euler's method
%   MaxStep         -- Estimates maximum step based on diffusion and force
%   TwoStepMilstein -- Estimates error using two half steps
%
% Dynamics Systems
%   Stokes          -- System without inertia
%   Newtonian       -- System with inertia

  % What do we do about walls, multiple particles?
  % We want adapters for beam+Tmatrix as well as explicit force equations,
  % and explicit optical potential (conservative).

%
% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

