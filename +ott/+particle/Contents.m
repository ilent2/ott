% ott.particle provides representations of particles.
%
% The Particle class aims to unify particle related properties such as
% drag and optical scattering methods.  At the moment this package
% only provides two Particle instances: Fixed and Variable.  Fixed
% stores a T-matrix, drag and shape instance and simplifies scattering
% and visualisation tasks.  Variable re-calculates the T-matrix and drag
% according to a user supplied method when the particle properties change.
%
% In a future version this interface may change to support other scattering
% methods, and/or more optimal T-matrix calculation methods (such as only
% building required columns of the T-matrix depending on the incident beam).
%
% Contents
%   Particle   - Base class for particles in optical tweezers simulations.
%   Fixed      - A particle with stored drag/tmatrix properties.
%   Variable   - A particle whose drag/tmatrix are automatically recomputed.
%
% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.
