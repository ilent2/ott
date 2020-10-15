% ott.drag Methods to calculate drag
%
% This package contains methods to calculate Stokes drag for spherical
% and arbitrary shaped particles in different configurations.
%
% Generic classes
%   Stokes         - Base class for Stokes drag tensors
%   StokesData     - Storage class for stokes drag tensors
%
% Free particle:
%   StokesSphere   - Drag tensor for a sphere with Stokes Drag
%   StokesCylinder - Drag tensor for a slender cylinder.
%   StokesLambNn   - Star shaped particles using pre-trained NN.
%   StokesLambPm   - Star shaped particles using Lamb series/point matching.
%
% Wall effects:
%   FaxenSphere    - Faxen's sphere corrections for movement near a plane.
%   PadeSphere     - Pade's sphere corrections for movement near a plane.
%   ChaouiSphere   - Choui's sphere corrections for movement near a plane.
%   EccentricSpheresNn - Eccentric sphere using Gibson's NN approach.
%
% Folders
%   private        - Functions used by drag calculation classes
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

