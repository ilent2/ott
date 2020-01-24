% ott.drag Methods to calculate drag
%
% Base class
%   Stokes              - Base class for Stokes drag tensors
%
% Free particle:
%   StokesSphere        - Drag tensor for a sphere with Stokes Drag
%   StokesLambNn        - Calculate stokes drag for star shaped particles using pre-trained NN.
%   StokesLambPm        - Calculate drag coefficients using Lamb series and point matching.
%
% Wall effects:
%   FaxenSphere         - Stokes drag with Faxen's corrections for movement near a plane.
%   EccentricSpheresNn  - Calculate drag on an eccentric sphere using Gibson's NN approach.
%
% Folders
%   private             - Functions used by drag calculation classes
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

