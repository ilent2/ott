% ott.bsc Beam shape coefficient classes (beams)
%
% Descriptions of beams in a vector spherical wave function basis.
% These classes encapsulate the VSWF coefficients.
%
% Classes
%   Bsc            - class representing beam shape coefficients
%   Scattered      - A Bsc instance describing scattered beams
%
%   Pointmatch     - base class for BSC generated using point matching
%   PmGauss        - provides HG, LG and IG beams using point matching method
%   PmParaxial     - calculate representation from farfield/paraxial beam
%   PmAnnular      - generate a beam with an annular far-field profile
%
%   Plane          - Plane wave
%
%   Bessel         - Bessel and Bessel-like beams with OAM
%   BesselAnnular  - generate a beam with an annular far-field profile

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

